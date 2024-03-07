%% Set up to run Jane task analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'S367';
condition = 'HL';
location = 2; 

if location == 1 % School
    prefix = 'F:\';
elseif location == 2 % SNAPlab
    prefix = 'E:\'; 
elseif location == 0 % Mac
    prefix = '/Volumes/SNH/';
end 

uname = 'samhauser';
suffix = ['THESIS/Pitch_Diagnostics_Data/Jane/Human/', condition,'/',subj];
datapath = [prefix,suffix];

Jane_Analysis;

