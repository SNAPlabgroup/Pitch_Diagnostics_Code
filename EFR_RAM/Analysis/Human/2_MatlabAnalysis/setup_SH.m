%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

uname = 'samhauser';

condition = 'YNH';
subj = 'S362';

prefix = ['/Users/samhauser/Library/CloudStorage/Box-Box/SNAPlab Data Archive/Pitch_Diagnostics_SH_AS/EFR_RAM/Human/'];
suffix = [condition,'/',subj,'/Preprocessed'];
datapath = [prefix,suffix];

processSubject;