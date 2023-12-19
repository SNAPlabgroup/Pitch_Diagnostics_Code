%Here's where you can define your own parameters for input/output
%directories.

clear;

condition = 'CA_2wksPost';
subj = 'Q431';
loc = 1; 

if loc == 1
    prefix = ['F:/', 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
    filesep, 'Chin', filesep]; 
elseif loc == 0
    prefix = ['/Volumes/SNH/','THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
    filesep, 'Chin', filesep]; 
end

suffix = [ condition,filesep,subj];
datapath = [prefix,suffix];


processChin;