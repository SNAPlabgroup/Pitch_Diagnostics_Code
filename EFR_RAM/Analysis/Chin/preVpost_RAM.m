

clear;

subj = 'Q443';
conditions = {'Baseline', 'Baseline'};

location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

for k = 1:length(conditions)
    condition = conditions{k};
    suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
        filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Processed'];
    datapath = [prefix,suffix];
    
    % Import Data
    cwd = pwd;
    cd(datapath)
    datafile = dir(fullfile(cd,[ subj, '_RAM_EFR_223_' condition, '*.mat']));
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
    
    t_all(k, :) = t; 
    T_env_all(k,:) = T_env; 
    PKS_all(k,:) = PKS;
    PLV_env_all(k,:) = PLV_env;
    f_all(k,:) = f;
    LOCS_all(k,:) = LOCS; 
    
    cd(cwd)
    
end

%% Plot:
blck = [0.25, 0.25, 0.25];
rd = [0.8500, 0.3250, 0.0980, 0.5];
figure;

%Spectral Domain
hold on;
title([subj,' | RAM - 25% Duty Cycle | ',condition],'FontSize',14);
plot(f_all(1,:),PLV_env_all(1,:),'Color',blck,'linewidth',1.5);
plot(LOCS_all(1,:),PKS_all(1,:),'*','Color',[.7, .7, .7],'MarkerSize',10,'LineWidth',2);

hold on;
title([subj,' | RAM - 25% Duty Cycle | ',condition],'FontSize',14);
plot(f_all(2,:),PLV_env_all(2,:),'Color',rd,'linewidth',1.5);
plot(LOCS_all(2,:),PKS_all(2,:),'*','Color',[.85, 0.32, 0],'MarkerSize',10,'LineWidth',2);


hold off;
ylim([0,1])
ylabel('PLV','FontWeight','bold')
xlabel('Frequency(Hz)','FontWeight','bold')

%Time Domain
xstart = .6;
xend = .9;
ystart = 0.6;
yend = .9;

axes('Position',[xstart ystart xend-xstart yend-ystart])
box on
hold on
plot(t_all(1,:), T_env_all(1,:),'Color',blck, 'LineWidth',2);
plot(t_all(2,:), T_env_all(2,:),'Color',rd, 'LineWidth',2);
xlim([0.3,.4]);
ylim([-2,2]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

set(gcf,'Position',[1557 538 560 420])