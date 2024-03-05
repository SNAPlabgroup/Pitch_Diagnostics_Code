%% Set up to run DPOAE analysis

%Here's where you can define your own parameters for input/output
%directories.


clear;

cwd = pwd; 
subj = 'Q443';
conditions = {'Baseline', 'PTS_2wksPost'};
location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

for k = 1:length(conditions)
    condition = conditions{k};
    suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'MEMR', ...
        filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Processed'];
    datapath = [prefix,suffix];
    
    % Import Data
 
    cd(datapath)
    datafile = dir(fullfile(cd,[subj, '_MEMR_WB_' condition, '*.mat']));
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
    
    power = mean(abs(data.res.MEM(:, data.res.ind)), 2);
    deltapow_all(k,:) = power - min(power);
    
    elic_all(k,:) = data.res.elicitor; 
  
    f_all(k,:) = data.res.freq;
    ind_all(k,:) = data.res.ind; 
    thresh_all(k,:) = data.res.threshold; 

    cd(pwd)
    
end

figure; 
hold on; 
%first
plot(elic_all(1,:), deltapow_all(1,:), 'ok-', 'linew', 2);

plot(elic_all(2,:), deltapow_all(2,:), 'or-', 'linew', 2);

%plot(elic_all(3,:), deltapow_all(3,:), 'ob-', 'linew', 2);

%labels
xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\Delta Absorbed Power (dB)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Pre', '2wksPost', 'location', 'Northwest')
%ymax = max(deltapow_all(:,:)+.05);
%ylim([0,ymax])
set(gca, 'XScale', 'log', 'FontSize', 14)
title(sprintf('WB-MEMR | %s', subj))
drawnow;

cd(cwd)