% Organize chin files
clear;

subj = 'Q427';
gender = 'M';
condition = 'PTS_2wksPost';
user = 'SH';
loc = 0;


%% Run
if strcmp(user, 'SH')
    if loc == 1
        prefix = ['F:\'];
    elseif loc == 0 % mac
        prefix = ['/Volumes/SNH/'];
    end
    suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, '_RawChinData', filesep, condition];
elseif strcmp(user, 'AS')
    if loc == 1
    end
end

alldatadir = [prefix suffix];

addpath('./convertData')
cwd = pwd;
cd(alldatadir)

chindirs = dir(fullfile(cd,['*',subj, '*']));

for dirNum = 1:length(chindirs)
    file = chindirs(dirNum).name;
    fprintf(sprintf('Loading directory: %s \n', file))
    cd(file)
    
    all_p_files = dir(fullfile(cd,['p*']));
    pfiles = {all_p_files.name};
    
    for p = 1:length(pfiles)
        
        if contains(pfiles(p), 'sweptDPOAE')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            load(char(pfiles(p)))
            [filename, data] = convertSweptDP(x, subj, file, condition,gender);
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition filesep subj filesep 'Preprocessed' filesep, filename];
            save(full_filename, 'data')
            %fprintf(sprintf('saved: %s\n', filename))
            clear x;
            
        elseif contains(pfiles(p), 'sweptSFOAE')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            load(char(pfiles(p)))
            [filename, data] = convertSweptSF(x, subj, file, condition,gender);
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition filesep subj filesep 'Preprocessed' filesep, filename];
            save(full_filename, 'data')
            %fprintf(sprintf('saved: %s\n', filename))
            clear x;
            
        elseif contains(pfiles(p), 'TEOAE')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'TEOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'TEOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'TEOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            
            load(char(pfiles(p)))
            [filename, data] = convertTEOAE(x, subj, file, condition,gender);
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'TEOAE', filesep, 'Chin', filesep, condition filesep subj filesep 'Preprocessed' filesep, filename];
            save(full_filename, 'data')
            %fprintf(sprintf('saved: %s\n', filename))
            clear x;
            
        elseif contains(pfiles(p), 'memr')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'MEMR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'MEMR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'MEMR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            load(char(pfiles(p)))
            [filename, data] = convertWBMEMR(x, subj, file, condition,gender);
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'MEMR', filesep, 'Chin', filesep, condition filesep subj filesep 'Preprocessed' filesep filename];
            save(full_filename, 'data')
            %fprintf(sprintf('saved: %s\n', filename))
            clear x;
            
        elseif contains(pfiles(p), '_dpoae')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'DPOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'DPOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'DPOAE', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('NEL DPOAE saved\n')
            
        elseif contains(pfiles(p), 'ABR') || contains(pfiles(p), 'calib_raw') || contains(pfiles(p), 'calib_inv')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('ABRs saved\n')
            
        elseif contains(pfiles(p), 'FFR_RAM')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('RAM saved\n')
            
        elseif contains(pfiles(p), 'PitchEFR') || contains(pfiles(p), 'f0103')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_p_files(p).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('Pitch EFRs saved\n')
        end
        
    end
    
    all_a_files = dir(fullfile(cd,['a*']));
    afiles = {all_a_files.name};
    
    for a = 1:length(afiles)
        
        if contains(afiles(a), 'ABR')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_a_files(a).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'ABR', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('ABR (a files) saved\n')
            
        elseif contains(afiles(a), 'FFR_RAM')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_a_files(a).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'EFR_RAM', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            % fprintf('RAM (a files) saved\n')
            
        elseif contains(afiles(a), 'PitchEFR') || contains(afiles(a), 'f0103')
            if ~exist(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file), 'dir')
                mkdir(strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                    filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            end
            copyfile(all_a_files(a).name, strcat(prefix, 'THESIS', filesep, 'Pitch_Diagnostics_Data', ...
                filesep, 'EFR_Pitch', filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Raw', filesep, file))
            %fprintf('PitchEFR (a files) saved\n')
        end
    end
    
    fprintf(sprintf('Saved files from directory: %s \n', file))
    cd(alldatadir)
end

cd(cwd)
