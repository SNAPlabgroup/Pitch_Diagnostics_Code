% Organize chin files

subj = 'Q412';
gender = 'M';
condition = 'Baseline';
user = 'SH';
loc = 0; 


%% Run
if strcmp(user, 'SH')
    if loc == 1
        prefix = ['F:\'];
    elseif loc == 0 % mac
        prefix = ['/Volumes/SNH/']; 
    end
end

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, '_RawChinData']; 

alldatadir = [prefix suffix]; 

addpath('./convertData')
cwd = pwd;
cd(alldatadir)

chindirs = dir(fullfile(cd,['*',subj, '*']));

for dirNum = 1:length(chindirs)
    file = chindirs(dirNum).name; 
    fprintf(sprintf('Loading directory: %s \n', file))
    cd(file)
    
    all_files = dir(fullfile(cd,['p*'])); 
    pfiles = {all_files.name}; 
    
    for p = 1:length(pfiles)
        
        if contains(pfiles(p), 'sweptDPOAE')
            load(char(pfiles(p)))
            [filename, data] = convertSweptDP(x, subj, file, condition,gender); 
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition filesep subj filesep filename];
            save(full_filename, 'data')
            fprintf(sprintf('saved: %s\n', filename))
            clear x; 
            
        elseif contains(pfiles(p), 'sweptSFOAE')
            load(char(pfiles(p)))
            [filename, data] = convertSweptSF(x, subj, file, condition,gender); 
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition filesep subj filesep filename];
            save(full_filename, 'data')
            fprintf(sprintf('saved: %s\n', filename))
            clear x; 
            
        elseif contains(pfiles(p), 'TEOAE')
            load(char(pfiles(p)))
            [filename, data] = convertTEOAE(x, subj, file, condition,gender); 
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'TEOAE', filesep, 'Chin', filesep, condition filesep subj filesep filename];
            save(full_filename, 'data')
            fprintf(sprintf('saved: %s\n', filename))
            clear x; 
            
        elseif contains(pfiles(p), 'memr')
            load(char(pfiles(p)))
            [filename, data] = convertWBMEMR(x, subj, file, condition,gender); 
            full_filename = [prefix 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'MEMR', filesep, 'Chin', filesep, condition filesep subj filesep filename];
            save(full_filename, 'data')
            fprintf(sprintf('saved: %s\n', filename))
            clear x; 
        end
        
    end
    cd(alldatadir)
end
cd(cwd)
