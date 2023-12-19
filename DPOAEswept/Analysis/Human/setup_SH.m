%% Run DPOAE Analysis Scripts for individual subject
% Human version

clear;

%% Enter information here: 
subj = 'S367';                    % e.g., 'Q440'
condition = 'HL';                 % e.g, 'HL', 'YNH', 'MANH'
location = 1;                     % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

%% Run Scripts
if location == 1 % School
    prefix = 'F:\';
elseif location == 2 % SNAPlab
    prefix = 'E:\'; 
elseif location == 0 % Mac
    prefix = '/Volumes/SNH/';
end 

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep,...
    'DPOAEswept', filesep, 'Human', filesep, condition, filesep, subj];

datapath = [prefix,suffix];

DPanalysis;

