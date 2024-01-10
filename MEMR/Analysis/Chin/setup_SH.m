%% Set up to run MEMR analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'Q441';
condition = 'Baseline';

user = 'SH'; 
loc = 0; 

uname = 'samhauser';

if strcmp(user, 'SH')
    if loc == 1
        prefix = ['F:\'];
    elseif loc == 0 % mac
        prefix = ['/Volumes/SNH/']; 
    end
end

suffix = ['THESIS/Pitch_Diagnostics_Data/MEMR/Chin/', condition,'/',subj, filesep, 'Preprocessed', filesep];
datapath = [prefix,suffix];

WBMEMR_Analysis;
