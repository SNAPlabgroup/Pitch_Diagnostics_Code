%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

condition = 'Baseline';
subj = 'Q421';

uname = 'asivapr';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/EFR_RAM/Chin/'];
suffix = [condition,'/',subj];

datapath = [prefix,suffix];

processChin;