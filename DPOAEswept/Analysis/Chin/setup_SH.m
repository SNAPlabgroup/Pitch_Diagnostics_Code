%% Set up to run SFOAE analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'Q426';
condition = 'Baseline';
location = 1; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end 

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'DPOAEswept', filesep, 'Chin', filesep, condition,filesep,subj];
datapath = [prefix,suffix];

suffix_calib = ['FPLcalib', filesep, 'Chin', filesep, condition, filesep, subj];
calibpath = [prefix, suffix_calib]; 

DPanalysis;
