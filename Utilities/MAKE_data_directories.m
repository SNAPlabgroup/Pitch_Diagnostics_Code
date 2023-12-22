%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make folders for new subject
% written by: S. Hauser
% August 1, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Enter info about new subject(s)
subj = 'Q426';
conditions = {'TTS_1dayPost', 'TTS_2wksPost'};
user = 'SH';
location = 0; % 0==laptop/mac; 1==work

%% Run
if strcmp(user, 'SH')
    if location == 0
        prefix = '/Volumes/SNH';
        fs = '/';
    elseif location == 1
        prefix = 'F:';
        fs = '\';
    elseif location == 2 % SNAPlab
        prefix = 'E:';
        fs = '\';
    else
        prefix = '/Volumes/';
    end
    datadir = [prefix fs 'THESIS' fs 'Pitch_Diagnostics_Data'];
    % else
    % Andrew directories here
end

measures_both = {'ABR', 'DPOAEswept', 'EFR_AM', 'EFR_Pitch', 'EFR_RAM', 'MEMR', 'SFOAEswept', 'TEOAE', 'FPLcalib'};
measures_chin = {'DPOAE'};
measures_human = {'ACC', 'AMdetection', 'ARDC', 'F0DL', 'Jane', 'MRT_tonecomplex', 'PTC'};

if subj(1) == 'Q'
    species = 'Chin';
else
    species = 'Human';
end
for s = 1:length(conditions)
    status = conditions{s}; 
    for m = 1:length(measures_both)
        meas = measures_both{m};
        if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj),'dir')) %Make subj folder and subfolders
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
            fprintf('Made %s dir for %s \n', meas, subj);
        end
        
        if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'),'dir'))
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
            fprintf('Made %s dir for %s (Preprocessed) \n', meas, subj);
        end
        if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'),'dir'))
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
            fprintf('Made %s dir for %s (Raw) \n', meas, subj);
        end
        if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'),'dir'))
            mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
            fprintf('Made %s dir for %s (Processed) \n', meas, subj);
        end
    end
    
    if strcmp(species,'Chin')
        for m1 = 1:length(measures_chin)
            meas = measures_chin{m1};
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
                fprintf('Made %s dir for %s \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
                fprintf('Made %s dir for %s (Preprocessed) \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
                fprintf('Made %s dir for %s (Raw) \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
                fprintf('Made %s dir for %s (Processed) \n', meas, subj);
            end
        end
    elseif strcmp(species,'Human')
        for m2 = 1:length(measures_human)
            meas = measures_human{m2};
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj),'dir')) %Make subj folder and all internal folders
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
                fprintf('Made %s dir for %s \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Preprocessed'));
                fprintf('Made %s dir for %s (Preprocessed) \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Raw'));
                fprintf('Made %s dir for %s (Raw) \n', meas, subj);
            end
            if (~exist(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'),'dir'))
                mkdir(strcat(datadir,fs,meas, fs, species, fs, status, fs, subj, fs, 'Processed'));
                fprintf('Made %s dir for %s (Processed) \n', meas, subj);
            end
        end
    end
    fprintf('All folders made for %s_%s \n', subj, status)
end



