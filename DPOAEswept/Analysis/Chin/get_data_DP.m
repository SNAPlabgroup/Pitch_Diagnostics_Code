function [oae, nf, f2, oaesum, centerfreq] = get_data_DP(prefix, subj, condition, AwakeFlag)


suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'DPOAEswept', ...
    filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Processed'];
datapath = [prefix,suffix];


% Import Data
cwd = pwd;
cd(datapath)

if ~AwakeFlag
    if exist('Sedated', 'dir')
        cd('Sedated')
    else
        % no sedated data
        fprintf(sprintf('No file for %s %s...Filling with NaNs!\n', subj, condition));
        oae=NaN(1,512);
        nf=NaN(1,512);
        f2=NaN(1,512);
        oaesum=NaN(1,9);
        centerfreq=NaN(1,9);
        cd(cwd);

        return
    end
end

datafile = dir(fullfile(cd,[subj, '_DPOAEswept_' condition, '*.mat']));
if length(datafile) < 1
    fprintf(sprintf('No file for %s %s...Filling with NaNs!\n', subj, condition));
    oae=NaN(1,512);
    nf=NaN(1,512);
    f2=NaN(1,512);
    oaesum=NaN(1,9);
    centerfreq=NaN(1,9);
elseif size(datafile,1) > 1
    checkDIR =uigetfile('.mat');
    load(checkDIR);
    file = checkDIR;
    load(file);
    oae = data.result.oae_full;
    nf = data.result.nf_full;
    f2= data.result.f2;
    oaesum = data.result.oae_summary;
    centerfreq = data.result.centerFreqs;
else
    load(datafile(1).name);
    file = datafile(1).name;
    load(file);
    oae = data.result.oae_full;
    nf = data.result.nf_full;
    f2= data.result.f2;
    oaesum = data.result.oae_summary;
    centerfreq = data.result.centerFreqs;
end



cd(cwd);


end