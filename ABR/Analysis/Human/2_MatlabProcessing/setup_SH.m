%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

uname = 'samhauser';

condition = 'HL';
subj = 'S367';

prefix = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/ABR/Human/'];
suffix = [condition,'/',subj,'/Preprocessed'];
datapath = [prefix,suffix];

processSubject;