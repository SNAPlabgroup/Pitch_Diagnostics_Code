%% Set up to run MEMR analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'S343';
condition = 'YNH';

uname = 'samhauser';
prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/MEMR/Human/'];
suffix = [condition,'/',subj];
datapath = [prefix,suffix];

WBMEMR_Analysis;
