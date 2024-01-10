%% Run DPOAE Analysis Scripts for individual subject
% Chinchilla version

clear;

%% Enter information here: 
subj = 'Q440';                  % e.g., 'Q440'
condition = 'Baseline';         % e.g., 'Baseline', 'CA_2wksPost'
location = 0;                   % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

%% Run Scripts
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

DPanalysis;
