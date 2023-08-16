%% Set up to run SFOAE analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'S360';
condition = 'YNH';
location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end 

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep,...
    'DPOAEswept', filesep, 'Human', filesep, condition,filesep,subj];
datapath = [prefix,suffix];

DPanalysis;

