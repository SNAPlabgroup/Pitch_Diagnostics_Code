%% Plotting the FPL INV 

subj = 'Q423';
condition = 'Baseline';         % e.g., 'Baseline', 'CA_2wksPost'
location = 0;                   % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';
% Run Scripts
if location == 1 % School
    prefix = 'F:\';
elseif location == 2 % SNAPlab
    prefix = 'E:\';
elseif location == 0 % Mac
    prefix = '/Volumes/SNH/';
end

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, ...
    'DPOAEswept', filesep, 'Chin', filesep, condition, filesep, subj];

datapath = [prefix,suffix];

% Import data
cwd = pwd;
cd([datapath, filesep, 'Preprocessed'])
datafile = dir(fullfile(cd,['sweptDPOAE_*.mat']));
if length(datafile) < 1
    fprintf('No file...Quitting!\n');
elseif size(datafile,1) > 1
    for j = 1:size(datafile,1)
        load(datafile(j).name);
        inv_f = data.FPL_inv.CalibData(:,1); 
        inv_level = data.FPL_inv.CalibData(:,2); 
        hold on;
        plot(inv_f, inv_level)
        
    end
else
    load(datafile(1).name);
    file = datafile(1).name;
end
legend(datafile.name, 'location', 'Northwest')
xlim([.5, 16])
set(gca, 'XScale', 'log')
cd(cwd);