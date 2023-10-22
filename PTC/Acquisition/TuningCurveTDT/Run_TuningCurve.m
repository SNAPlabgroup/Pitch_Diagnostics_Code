%% Tuning curve through Bekesy tracking of simultaneous masking thresholds.
%
% Pure tone pulses at fixed level are embeded in narrow-band noise that
% slowly sweeps either up from a low center frequency or down from a high
% center frequency.
%
% If spacebar is held down, noise level increases. Otherwise it decreases.
%
% The procedure and parameters are based on Sek et al., 2005.
%
% ----------------------------------------
% Copyright 2022- Hari Bharadwaj. All rights reserved.
% hari.bharadwaj@pitt.edu
%
% References:
% Sek, A., Alcantara, J., Moore, B. C., Kluk, K., & Wicher, A. (2005).
% Development of a fast method for determining psychophysical tuning
% curves. International Journal of Audiology, 44(7), 408-420.
%
% ----------------------------------------

try
    subj = input('Please subject ID:', 's');
    % Make directory to save results
    paraDir = './Results/';
    if(~exist(strcat(paraDir,'/',subj),'dir'))
        mkdir(strcat(paraDir,'/',subj));
    end
    respDir = strcat(paraDir,'/',subj,'/');
    
    %From 3AFC threshold data
    thresh = input('Enter threshold at 4 kHz in dB SPL:'); 

    fprintf(1, 'Press and hold SPACE bar if you DONT hear the beeps!\n');

    
    % Set some parameters about range of intensities and step sizes
    Lmax = 90;  % dB SPL
    Lmin = -20; % dB SPL
    Lstart = 50;  % dB SPL
    Lcurr = Lstart;
    Lrate = 2;  % dB per second (going up or down depending on key state)
    
    Ltone = thresh + 20;
    Lref = 114; % dB SPL expected from full-scale (-1 to 1) tone for ER-2
    
    
    %% Some setting up for keyboard methods.
    KbName('UnifyKeyNames'); % Just in case
    plevel = 2;
    
    oldlevel = Priority(plevel);
    escapeKey = KbName('ESCAPE');
    spaceKey = KbName('space');
    ListenChar(2);  % Don't show keyboard presses on console
    
    %% Loading stimuli that will be played and setting frame size
    
    % Contains the following variables
    % fs        : Sampling rate (Hz)
    % freqseq   : Local noise center frequency around each sample (Hz)
    % noise     : Actual noise waveform (full scale)
    % tone      : Pulsed tone waveform (full scale)
    % bw        : Local bandwidth of the noise (Hz)
    % startF    : Starting frequency of noise (Hz)
    % endF      : Ending frequency of noise (Hz)
    load('tuningStims.mat');
    
    
    frameDur_samps = 2048;
    
    noise = noise * sqrt(2) / rms(noise);
    % Lref is based on full scale tone which has rms of 1/sqrt(2)
    % Adjust noise accordingly

    
    Nframes = floor(numel(noise)/frameDur_samps);
    
    %% Create variables to hold tracking results
    % The code writes down the local frequency and noise level in each
    % frame. Reversal points can be extracted offline as needed.
    results.ftrack = zeros(Nframes, 1);
    results.Ltrack = zeros(Nframes, 1);
    
    % Write other useful parameters to results sstructure
    results.bw = bw;
    results.fs = fs;
    results.startF = startF;
    results.endF = endF;
    results.Lmax = Lmax;
    results.Lmin = Lmin;
    results.Ltone = Ltone;
    results.Lrate = Lrate;
    results.Lstart = Lstart;
    
    %% Calculate some level control parameters
    
    frameDur = frameDur_samps / fs; % In seconds
    dB_per_frame = frameDur * Lrate; % How much will level change is one frame
    dB_upward = linspace(0, dB_per_frame, frameDur_samps);
    dB_downward = linspace(0, -dB_per_frame, frameDur_samps);
    
    
    %% Initialize TDT RZ6 for blocked (frame-by-frame) use.
    fig_num=99;
    USB_ch=1;
    FS_tag = 3;
    Fs = 48828.125;
    [f1RZ,RZ,~]=load_play_circuit_block(FS_tag,fig_num,USB_ch);
    invoke(RZ, 'SetTagVal', 'blocksize',frameDur_samps);
    frameDur = frameDur_samps / Fs;

    dropL = 120;
    dropR = 0;
    invoke(RZ, 'SetTagVal', 'attA', dropL);
    invoke(RZ, 'SetTagVal', 'attB', dropR);
    invoke(RZ, 'SetTagVal', 'chanselect', 2);
    
    %% Start keeping time and playing audio
    t0 = tic;

    % Load first frame and start playout
    idx =  1:frameDur_samps;
    buff = noise(idx) .* db2mag(Lcurr - Lref) + ...
        tone(idx) * db2mag(Ltone - Lref);
    invoke(RZ, 'WriteTagVEX', 'datain', 0, 'F32', buff);
    invoke(RZ, 'SoftTrg', 1);

    for k = 2:Nframes
        prev = (k - 2) *frameDur_samps;
        start = (k - 1)*frameDur_samps;
        idx =  start + (1:frameDur_samps);
        
        % Calculate level-specific multiplication factors
        if Lcurr < Lmax
            upscale_factor = db2mag(Lcurr - Lref + dB_upward);
        else
            upscale_factor = db2mag(Lcurr - Lref);
        end
        
        if Lcurr > Lmin
            downscale_factor = db2mag(Lcurr - Lref + dB_downward);
        else
            downscale_factor = db2mag(Lcurr - Lref);
        end
        
        
        % Increase or decrease level depending on key state
        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
        keyCode = find(keyCode, 1);
        if keyIsDown
            if keyCode == spaceKey
                buff = noise(idx) .* downscale_factor + ...
                    tone(idx) * db2mag(Ltone - Lref);
                Lcurr = Lcurr - dB_per_frame;
                results.ftrack(k) = mean(freqseq(idx));
                results.Ltrack(k) = Lcurr + dB_per_frame / 2;
            end
            if keyCode == escapeKey
                break;
            end
        else
            buff = noise(idx) .* upscale_factor + ...
                tone(idx) * db2mag(Ltone - Lref);
            Lcurr = Lcurr + dB_per_frame;
            results.ftrack(k) = mean(freqseq(idx));
            results.Ltrack(k) = Lcurr - dB_per_frame / 2;
        end

        invoke(RZ, 'WriteTagVEX', 'datain', start, 'F32', buff);


        buffpos = invoke(RZ, 'GetTagVal', 'buffpos');
        while buffpos < prev
            WaitSecs(0.1*frameDur);
            buffpos = invoke(RZ, 'GetTagVal', 'buffpos');
        end

        Lcurr = max(Lcurr, Lmin);
        Lcurr = min(Lcurr, Lmax);
    end
    toc(t0);

    
    datetag = char(string(datetime('now')));
    datetag(strfind(datetag,' ')) = '_';
    datetag(strfind(datetag,':')) = '_';
    saveName = [respDir, '/', subj, '_Tuning_', datetag, '.mat'];
    save(saveName, 'results');

    %Clearing I/O memory buffers:
    invoke(RZ,'ZeroTag','datain');
    close_play_circuit(f1RZ,RZ);
    
    ListenChar(0);
    Priority(oldlevel);
catch me
    %Clearing I/O memory buffers:
    invoke(RZ,'ZeroTag','datain');
    close_play_circuit(f1RZ,RZ);
    ListenChar(0);
    Priority(oldlevel);
    rethrow(me);
end
