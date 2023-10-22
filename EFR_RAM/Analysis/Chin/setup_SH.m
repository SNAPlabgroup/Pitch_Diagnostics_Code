%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

condition = 'Baseline';
subj = 'Q422';


prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/EFR_RAM/Chin/'];
suffix = [condition,'/',subj];
datapath = [prefix,suffix];


processChin;