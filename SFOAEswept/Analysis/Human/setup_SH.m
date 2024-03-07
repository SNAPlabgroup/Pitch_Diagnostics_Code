%% Set up to run SFOAE analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'S368';
condition = 'HL';

location = 0; 

if location == 1 % School
    prefix = 'F:\';
elseif location == 2 % SNAPlab
    prefix = 'E:\'; 
elseif location == 0 % Mac
    prefix = '/Volumes/SNH/';
end 

uname = 'samhauser';

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep,...
    'SFOAEswept', filesep, 'Human', filesep, condition,filesep,subj];
datapath = [prefix,suffix];


SFanalysis;

