%Setup file for directory information (ABR threshold)

close all;
clear;

condition = 'Baseline';
subj = 'Q428';

uname = 'asivapr';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/ABR/Chin/'];
suffix = [condition,'/',subj];
datapath = [prefix,suffix];

freqs = [500,1e3,2e3,4e3,8e3];
% freqs = [8000];
ABR_audiogram_chin;

