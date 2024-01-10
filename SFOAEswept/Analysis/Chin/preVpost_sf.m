

clear;

subj = 'Q443';
conditions = {'Baseline', 'Baseline', 'Baseline','Baseline','Baseline'};
location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

for k = 1:length(conditions)
    condition = conditions{k};
    suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, ...
        'SFOAEswept', filesep, 'Chin', filesep, condition, filesep, subj, ...
        filesep, 'Processed'];
    datapath = [prefix,suffix];
    
     AwakeFlag = 1; 
    
    % Import Data
    cwd = pwd;
    cd(datapath)
    if exist('Sedated', 'dir')
        sed = questdlg('Awake or Sedated?', 'Sedated?', 'Awake', 'Sedated', 'Awake'); 
        if strcmp(sed, 'Sedated')
            cd('Sedated')
            AwakeFlag = 0; 
        end  
    end


    datafile = dir(fullfile(cd,[subj, '_SFOAEswept_' condition, '*.mat']));
    if length(datafile) < 1
        fprintf('No file...Quitting!\n');
    elseif size(datafile,1) > 1
        checkDIR =uigetfile('.mat');
        load(checkDIR);
        file = checkDIR;
    else
        load(datafile(1).name);
        file = datafile(1).name;
    end
    
    load(file);
    
    oae_all(k, :) = data.result.oae_full;
    nf_all(k,:) = data.result.nf_full;
    f_all(k,:) = data.result.f;
    oaesum_all(k,:) = data.result.oae_summary;
    centerfreq_all(k,:) = data.result.centerFreqs;
    
    cd(cwd)
    
end

figure;
hold on;


colors = {'k', 'r', 'b', 'g', 'y'}; 

for k = 1:length(conditions)
    semilogx(f_all(k,:), oae_all(k,:), 'Color', colors{k}, 'linew', 2)
    %legend(conditions)
end

for k=1:length(conditions)
    semilogx(f_all(k,:), nf_all(k,:), '--', 'linew', 1.5, 'Color', colors{k})
    semilogx(centerfreq_all(k,:), oaesum_all(k,:), 'o', 'linew', 2, 'MarkerSize', 8, 'MarkerFaceColor', colors{k}, 'MarkerEdgeColor', colors{k})
end


%legend('pre-oae', 'pre-nf', '', '1d post-oae', '1d post-nf', '', '2w post-oae', '2w post-nf', '','location', 'northwest')
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 60])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('F Frequency (kHz)', 'FontWeight', 'bold')
title(['SFOAE | ',subj], 'FontSize', 16);

cd(cwd)
