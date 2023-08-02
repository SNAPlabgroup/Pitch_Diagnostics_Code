%% Set up to run SFOAE analysis

%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

subj = 'S361';
condition = 'YNH';

uname = 'samhauser';
prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/DPOAEswept/Human/'];
suffix = [condition,'/',subj];
datapath = [prefix,suffix];

DPanalysis;

