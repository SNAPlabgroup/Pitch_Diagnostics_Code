%% Set up to run TEOAE analysis

%Here's where you can define your own parameters for input/output
%directories.


clear;

subj = 'Q431';
condition = 'CA_2wksPost';
location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end 

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'TEOAE', filesep, 'Chin', filesep, condition, filesep, subj];
datapath = [prefix,suffix]; 

TEOAE_Analysis;
