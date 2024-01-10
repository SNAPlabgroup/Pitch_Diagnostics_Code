%% Get Data from SNAP lab computer

subj = 'S367';
condition = 'HL';
user = 'SH';

%% Run
if strcmp(user, 'SH')
    fs = '\';
    datadir = ['E:' fs 'THESIS' fs 'Pitch_Diagnostics_Data'];
    % else
    % Andrew directories here
end

% SNAP dir
snapdir = ['C:\Experiments\'];
cdir = pwd;

measdirs = {'Sam\DPOAEswept', 'Sam\SFOAEswept', 'Sam\Threshold_ForTuning', 'Sam\TuningCurveTDT' ...
    'Sam\TEOAE', 'Sam\MRT_tonecomplex', 'Sam\Jane', 'Sam\FrequencyTuning_SekMoore', ...
    'Sam\AM', 'Sam\WBMEMR', 'Andrew\Threshold_ForTuning', 'Andrew\FrequencyTuning_SekMoore', ...
    'Andrew\spectL_F0DL'};
measfolder = {'DPOAEswept', 'SFOAEswept', 'PTC', 'PTC',...
    'TEOAE', 'MRT_tonecomplex', 'Jane', 'PTC', ...
    'AMdetection', 'MEMR', 'PTC', 'PTC', ...
    'F0DL'};
for i = 1:length(measdirs)
    measdir = [snapdir measdirs{i}];
    meas = measfolder{i};
    if exist(strcat(measdir, '\Results\', subj),'dir')
        cd(strcat(measdir, '\Results\', subj))
        files = dir('*.mat');
        for f = 1:length(files)
            copyfile(files(f).name, strcat(datadir, fs, meas, fs, 'Human', fs, condition, fs, subj, fs, 'Raw'))
            fprintf('Copy file: %s\n', files(f).name)
        end
        cd(snapdir)
    end
end

if exist(strcat(snapdir, '\FPLclick\EARCAL\', subj),'dir')
    cd(strcat(snapdir, '\FPLclick\EARCAL\', subj));
    fpl_files = dir('*.mat');
    for f = 1:length(fpl_files)
        copyfile(fpl_files(f).name,  strcat(datadir, fs, 'FPLcalib', fs, 'Human', fs, condition, fs, subj, fs, 'Raw'))
        fprintf('Copy file: %s\n', fpl_files(f).name)
    end
    cd(snapdir)
end

cd(cdir)