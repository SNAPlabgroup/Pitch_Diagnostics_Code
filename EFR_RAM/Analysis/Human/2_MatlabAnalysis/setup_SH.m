%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

uname = 'samhauser';

condition = 'HL';
subj = 'S366';

prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/EFR_RAM/Human/'];
suffix = [condition,'/',subj,'/Preprocessed'];
datapath = [prefix,suffix];

processSubject;