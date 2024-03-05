function [power, deltapows, elicitors, thresh] = get_data_memr(prefix, subj, condition)


suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'MEMR', ...
    filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Processed'];
datapath = [prefix,suffix];


% Import Data
cwd = pwd;
cd(datapath)
datafile = dir(fullfile(cd,[subj, '_MEMR_WB_' condition, '*.mat']));

if length(datafile) < 1
    fprintf(sprintf('No file for %s %s...Filling with NaNs!\n', subj, condition));
    
    power = nan(1,11); 
    deltapows = nan(1,11);
    elicitors = nan(1,11);
    thresh = nan(1,1); 
    
elseif size(datafile,1) > 1
    checkDIR =uigetfile('.mat');
    load(checkDIR);
    file = checkDIR;
    load(file);
    power = mean(abs(data.res.MEM(:, data.res.ind)), 2)';
    deltapows = power - min(power);
    elicitors = data.res.elicitor';
    thresh = data.res.threshold; 
else
    load(datafile(1).name);
    file = datafile(1).name;
    load(file);
    power = mean(abs(data.res.MEM(:, data.res.ind)), 2);
    deltapows = power - min(power);
    elicitors = data.res.elicitor;
    thresh = data.res.threshold; 
end



cd(cwd);


end
