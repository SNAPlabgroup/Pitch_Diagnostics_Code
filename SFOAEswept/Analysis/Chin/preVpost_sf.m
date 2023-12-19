

clear;

subj = 'Q430';
conditions = {'Baseline', 'CA_2wksPost'};
location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

for k = 1:length(conditions)
    condition = conditions{k};
    suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'SFOAEswept', filesep, 'Chin', filesep, condition, filesep, subj];
    datapath = [prefix,suffix];
    
    % Import Data
    cwd = pwd;
    cd(datapath)
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

% for p = 1:size(oaesum_all,1)
% semilogx(f_all(p,:), oae_all(p,:), 'linew', 2)
% semilogx(f_all(p,:), nf_all(p,:), '--', 'linew', 1.5)
% semilogx(centerfreq_all(p,:), oaesum_all(p,:), 'o', 'linew', 2, 'MarkerSize', 8)
% end

semilogx(f_all(1,:), oae_all(1,:), 'Color', 'k', 'linew', 2)
semilogx(f_all(1,:), nf_all(1,:), '--', 'linew', 1.5, 'Color', 'k')
semilogx(centerfreq_all(1,:), oaesum_all(1,:), 'o', 'linew', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')

semilogx(f_all(2,:), oae_all(2,:), 'Color', 'r', 'linew', 2)
semilogx(f_all(2,:), nf_all(2,:), '--', 'linew', 1.5, 'Color', 'r')
semilogx(centerfreq_all(2,:), oaesum_all(2,:), 'o', 'linew', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')

semilogx(f_all(3,:), oae_all(3,:), 'Color', 'b', 'linew', 2)
semilogx(f_all(3,:), nf_all(3,:), '--', 'linew', 1.5, 'Color', 'b')
semilogx(centerfreq_all(3,:), oaesum_all(3,:), 'o', 'linew', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')

semilogx(f_all(4,:), oae_all(4,:), 'Color', 'g', 'linew', 2)
semilogx(f_all(4,:), nf_all(4,:), '--', 'linew', 1.5, 'Color', 'g')
semilogx(centerfreq_all(4,:), oaesum_all(4,:), 'o', 'linew', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g')


%legend('pre-oae', 'pre-nf', '', '1d post-oae', '1d post-nf', '', '2w post-oae', '2w post-nf', '','location', 'northwest')
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 60])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('F Frequency (kHz)', 'FontWeight', 'bold')
title(['SFOAE | ',subj], 'FontSize', 16);

cd(cwd)
