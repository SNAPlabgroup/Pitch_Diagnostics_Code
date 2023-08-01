clear all; close all hidden; clc; %#ok<CLALL>

% NOTE: db_main and db_flank are re-scaled to make total tone complex level
% 80 db SPL
% ranks = [4,7,10,13,16]; %from Oxenham paper
% ranks = [4,6,8,10,12,14];
ranks = 2:2:12;
ntrials = 2500; % This is the number of trials per rank
right = 1;
F0 = 103;
ALT_phase = 1;
nharms = 4;
db_main = 70;
db_flank = 64;
level = 80; %tone complex db SPL
noise_on = 0;
rove = 0;
stim_dur = .2;
tone_rms = 0.020;
ramp = 0.02;
phi = ones(nharms+2, 1) * pi/2; % Sin for all harmonics
rove = 0; %No need to rove here, stick with normal harmonic rank

if ALT_phase
    phi(1:2:end) = 0;
    phi = phi + pi/4; %add to make sure starting tone complex at 0
end

try
    % Loading random generator seed and state so that anything generated
    % randomly will be regenerated with the same realization everytime
    load('s.mat');
    rng(s);
    
    % Initialize play circuit
    fig_num=99;
    USB_ch=1;
    IAC = -1;
    FS_tag = 3;
    Fs = 48828.125;
    
    [f1RZ,RZ,FS]=load_play_circuit(FS_tag,fig_num,USB_ch,0,IAC);
    
    % Experiment parameters
    tot_trials = ntrials*length(ranks);
   
    isi = 0.18; % Average interstimulus interval (only approx guarenteed)

    % Some cushion for time taken to write each trial to buffer
    % Stim has to be short enough to be written in this time, typically 50
    % ms should be enough for most stims.
    bufferwritetime = 0.02;
    jitter = 0.05; % Maximum value of uniformly distributed jitter in ISI
    
    % Send trigger to EEG to start saving
    invoke(RZ, 'SetTagVal', 'trigval',253);
    invoke(RZ, 'SoftTrg', 6);
    
    % Wait at least 3 seconds before presenting any stims (e.g., for
    % filter transients, etc.)
    WaitSecs(3.0);
    
    % Load or generate a stimulus with the intent of playing it in both
    % polarities. Replace with other stims here.
    % Must be mono (i.e., size Nsamples-by-1) with a variable called 'y'.
    % This can also be generated within the main loop if each trial has a
    % different stimulus, but depending on stim generation code, in some
    % cases that might slowdown the overall timing (but not affect the
    % synchrony of triggers).
    % The sampling rate has to Fs = 48828.125 Hz based on earlier hardcoded
    % properties of this template script.
    
    rlist = repmat(1:length(ranks),1,tot_trials); 
    rlist = ranks(rlist); %rlist is a vector of tot_trials harmonic ranks
    rlist = rlist(randperm(numel(rlist))); %randomize elements
    
    % Using jitter to make sure that background noise averages out across
    % trials. We use jitter that is random between 0 and 'jitter' seconds.
    % Average duration added by the jitter is jitter/2 seconds
    jitlist = rand(tot_trials, 1)*jitter;
    
    if isi < 0.05
        error('Interstimulus interval too short, cannot continue');
    end
    
    % Keep track of time if needed
    tstart = tic;
    
    tonecount = 0;
    for j = 1:tot_trials
        
        curr_rank = rlist(j);
        
%         [x, sigrms] = makeComplexTone_Mehta_F0DL(F0, stim_dur, Fs, db_main, db_flank, curr_rank, nharms, ramp, phi, noise_on, rove);
        [x, sigrms] = make_F0DL_stim(F0, stim_dur, Fs, db_main, db_flank, tone_rms, curr_rank, nharms, ramp, phi, noise_on, rove);
        
        tonecount = tonecount + 1;

            p = (rand < 0.5); % 0 for plus, 1 for minus
            
            if(p==0)
                y = x;
                poloffset = 0;
            else
                y = x * -1;
                poloffset = 100;
            end
            
            stimrms = sigrms; %signal, no noise

            chanL = y;
            chanR = y;
            
            stimTrigger = curr_rank + poloffset;
            
            jit = jitlist(tonecount);
            
            stimlength = numel(y); % Recalculate, just in case
            
            %---------------------------------------------------------
            % Stimulus calibration calculations based on known
            % hardware and known .rcx circuit properties
            %---------------------------------------------------------
            % ER-2s give 100dB SPL for a 1kHz tone with a 1V-rms drive.
            % BasicPlay.rcx converts +/- 1 in MATLAB to +/- 5V at TDT
            % output, so you get 114 dB SPL for MATLAB rms of 1.
            %
            % So if we want a signal with rms of "stimrms" to be X dB,
            % then you have to attenuate in harwdware by below amount:
            
            % If switching to monaural, drop the unused ear by 120 dB flat
            % Also change overall stim level (because no loudness
            % summation)
            dropL = 114 - level + db(stimrms);
            dropR = 114 - level + db(stimrms); % Assumes diotic
            
            %to be absolutely sure monaural
            if right
                dropL = 120;
            else
                dropR = 120;
            end
            
            invoke(RZ, 'SetTagVal', 'trigval', stimTrigger);
            invoke(RZ, 'SetTagVal', 'nsamps', stimlength);
            %write to buffer left ear
            invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', chanL);
            %write to buffer right ear
            invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', chanR);
            
            %setting analog attenuation L
            invoke(RZ, 'SetTagVal', 'attA', dropL);
            %setting analog attenuation R
            invoke(RZ, 'SetTagVal', 'attB', dropR);
            
            % Just giving time for data to be written into buffer
            WaitSecs(bufferwritetime);
            
            %Start playing from the buffer:
            invoke(RZ, 'SoftTrg', 1); %Playback trigger
            fprintf(1,' Sequence %d/%d, tonecount = %d\n', j, tot_trials,...
                tonecount);
    
            WaitSecs(stim_dur + isi + jit);
        
    end
    invoke(RZ,'ZeroTag','datainL');
    invoke(RZ,'ZeroTag','datainR');
    
    toc(tstart); % Just to help get a sense of how long things really take
    
    %Clearing I/O memory buffers:

    % Have at least 3 seconds of no stim EEG data at the end
    WaitSecs(3.0);
    
    % Send trigger to EEG computer asking it to stop saving
    invoke(RZ, 'SetTagVal', 'trigval', 254);
    invoke(RZ, 'SoftTrg', 6);
    
    close_play_circuit(f1RZ,RZ);
    fprintf(1,'\n Done with data collection!\n');
    
catch me
    close_play_circuit(f1RZ,RZ);
    rethrow(me);
end
