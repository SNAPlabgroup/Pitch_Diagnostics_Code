%% CLICK OAE using traditional windowing method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Samantha Hauser, AuD
% Modified from: Hari Bharadwaj, PhD (SNAP Lab)
% Created: November 2021
% Last revision: 16-Sep-2023 (added metadata saving)
%
% References:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up data storage and subject info

% Measure-General info
info.measure = 'TEOAE';
info.version = 'v01';

% Visit info
if exist('C:\Experiments\Sam\current_visit.mat','file')
    load('C:\Experiments\Sam\current_visit.mat', 'visit')
    ask = questdlg(sprintf('Is this subject %s?', visit.subj.ID), ...
        'Check Subject', 'Yes', 'No', 'No');
else
    ask = 'No';
end

if strcmp(ask, 'No')
    cd ..
    startVisit
    cd(info.measure)
end

subj = visit.subj;
info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;

% Get ear info
subj.ear = questdlg('Which ear?', 'Ear', 'L', 'R', 'R');

% Get date/time
datetag = datestr(clock);
info.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '-';

% Make directory to save results
paraDir = 'C:\Experiments\Sam\TEOAE\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

tic; 
try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT and specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % load stimulus parameters
    click = clickSetDefaults();

    driver = 1; 
    drivername = strcat('Ph',num2str(driver));
    click.drivername = drivername;
    
    % Make click
    vo = clickStimulus(click.BufferSize + click.StimWin);
    buffdata = zeros(2, numel(vo));
    buffdata(driver, :) = vo; % The other source plays nothing
    click.vo = vo;
    odd = 1:2:click.Averages;
    even = 2:2:click.Averages;

    drop = click.Attenuation;
    dropOther = 120;
    fprintf('Starting Stimulation...\n')
    
    if driver == 1
        vins = playCapture2(buffdata, card, click.Averages, ...
            click.ThrowAway, drop, dropOther, 1);
    else
        vins = playCapture2(buffdata, card, click.Averages, ...
            click.ThrowAway, dropOther, drop, 1);
    end
    
    %compute the average
    vins = vins(:, (click.StimWin+1):(click.StimWin + click.RespDur)); % Remove stimulus by windowing
    
    if click.doFilt
        % High pass at 200 Hz using IIR filter
        [b, a] = butter(4, 200 * 2 * 1e-3/click.SamplingRate, 'high');
        vins = filtfilt(b, a, vins')';
    end
    
    vavg_odd = trimmean(vins(odd, :), 20, 1);
    vavg_even = trimmean(vins(even, :), 20, 1);
    rampdur = 0.2e-3; %seconds
    Fs = click.SamplingRate/2 * 1e3;
    resp.vavg = rampsound((vavg_odd + vavg_even)/2, Fs, rampdur);
    resp.noisefloor = rampsound((vavg_odd - vavg_even)/2, Fs, rampdur);
    
    Vavg = rfft(resp.vavg);
    Vavg_nf = rfft(resp.noisefloor);
    
    % Apply calibrations to convert voltage to pressure
    % For ER-10X, this is approximate
    mic_sens = 50e-3; % mV/Pa. TO DO: change after calibration
    mic_gain = db2mag(gain + 6); % +6 for balanced cable
    P_ref = 20e-6;
    DR_onesided = 1;
    factors = DR_onesided * mic_gain * mic_sens * P_ref;
    resp.output_Pa_per_20uPa = Vavg / factors; % unit: 20 uPa / Vpeak
    resp.noise_Pa_per_20uPa = Vavg_nf / factors;
    
    resp.freq = 1000*linspace(0,click.SamplingRate/2,length(Vavg))';
    
    
    %% Plot data
    PlotResults_simple;
    
    %% Save Ear Measurements
    data.info = info;
    data.stim = click;
    data.resp = resp; 
    data.info.subj = subj;
    data.resp.allTrials = vins;
    data.resp.testDur_s = toc;
    
    save(fname,'data');
    
    %% Close TDT, ER-10X connections etc. and cleanup
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    
catch me
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    rethrow(me);
end