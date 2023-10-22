%% Set up to run TEOAE analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'Q428';
condition = 'Baseline';

uname = 'samhauser';
prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/TEOAE/Chin/'];
suffix = [condition,'/',subj];
datapath = [prefix,suffix];

TEOAE_Analysis;
