%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

condition = 'Baseline';
subj = 'Q425';
loc = 1; 

if loc == 1
    prefix = ['F:/']; 
elseif loc == 2
    prefix = ['/Volumes/SNH/']; 
end

suffix = ['THESIS/Pitch_Diagnostics_Data/EFR_RAM/Chin/' condition,'/',subj];
datapath = [prefix,suffix];


processChin;