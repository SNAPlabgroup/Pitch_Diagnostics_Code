clear all; close all hidden; clc; %#ok<CLALL>

fig_num=99;
USB_ch=1;
IAC = -1;
FS_tag = 3;
Fs = 48828.125;

[f1RZ,RZ,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,0,IAC);

% Stimulus Parameters
dur = 2.0;
stimpe = 0.95; % Amplitude of each click
clickdur = 103; % microsecond (5 samples at 48.828125 Hz)
clickrate = 21; %Hz
jitwin = 0.04; 
levels = [93, 115];

% Experiment parameters
nconds = numel(levels);
ntrials = 240; % About 40 mins
order = repmat(1:nconds, 1, ntrials);
order = order(randperm(nconds*ntrials));
triglistL = [6, 9, 10, 12];
triglistR = [96, 144, 160, 192];

% Pause Off
invoke(RZ, 'SetTagVal', 'trigvalL',253);
invoke(RZ, 'SetTagVal', 'onsetdel',100);
invoke(RZ, 'SoftTrg', 6);

pause(2.0);



jitlist = rand(1, numel(order))*0.1;

t = 0:(1/Fs):(dur - 1/Fs);
tstart = tic;
for j = 1:numel(order)
    level = levels(order(j));
    for polarity = [1, -1]
        y = makeCABRstim_LR(clickdur, clickrate, dur, Fs, jitwin, polarity);

        chanL = y(1, :);
        chanR = y(2, :);
        

        
        jit = jitlist(j);
        
        stimlength = size(y, 2);
        stimTrigger = order(j);
        
        
        if polarity == -1
            stimTrigger = stimTrigger + nconds;
        end
        Ltrigger = triglistL(stimTrigger);
        Rtrigger = triglistR(stimTrigger);
        
        %-----------------
        % Why 111 for ER-2?
        %-----------------
        % ER-2s give about 100dB SPL for a 1kHz tone with a 1V-rms drive.
        % Max output is +/-5V peak i.e 3.54V-rms which is 11 dB higher.
        % Thus 111 dB-SPL is the sound level for tones when they occupy full
        % range.
        
        % Full range in MATLAB for a pure tone is +/- 1 which is 0.707 rms and
        % that corresponds to 111 dB SPL at the end. So if we want a signal
        % with rms sigrms to be x dB, then (111 - x) should be
        % db(sigrms/0.707).
        
        drop = 111 - level + 3 + db(stimpe); % The extra 3 for the 0.707 factor
        
        invoke(RZ, 'SetTagVal', 'trigvalL',Ltrigger);
        invoke(RZ, 'SetTagVal', 'trigvalR',Rtrigger);
        invoke(RZ, 'SetTagVal', 'nsamps', stimlength);
        invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', chanL); %write to buffer left ear
        invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', chanR); %write to buffer right ear
        
        invoke(RZ, 'SetTagVal', 'attA', drop); %setting analog attenuation L
        invoke(RZ, 'SetTagVal', 'attB', drop); %setting analog attenuation R
        
        WaitSecs(0.1); % Just giving time for data to be written into buffer
        %Start playing from the buffer:
        invoke(RZ, 'SoftTrg', 1); %Playback trigger
        fprintf(1,' Trial Number %d/%d, Level = %f peSPL, polarity = %d \n', j, ...
            numel(order), level, polarity);
        pause(dur + 0.1 + jit);
    end
end
toc(tstart);
%Clearing I/O memory buffers:
invoke(RZ,'ZeroTag','datainL');
invoke(RZ,'ZeroTag','datainR');
pause(3.0);

% Pause On
invoke(RZ, 'SetTagVal', 'trigvalL', 254);
invoke(RZ, 'SetTagVal', 'onsetdel',100);
invoke(RZ, 'SoftTrg', 6);

close_play_circuit(f1RZ,RZ);
fprintf(1,'\n Done with data collection!\n');

