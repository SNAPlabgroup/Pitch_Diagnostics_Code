%Author: Andrew Sivaprakasam
%Updated: July, 2023
%Purpose: Script to import/plot/apply additional processing to RAM_EFR
%files

%TODO: pass time bounds thru

close all
clear;

condition = 'YNH';
subj = 'S353';
fmod = 103;
harmonics = 4;

t_win = [.2,.9]; %signal window, ignoring onset/offset effects
filts = [60,4000];
%% Handles my Local File Structure/EXT drive

uname = 'sivaprakasaman';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/EFR_RAM/Human/'];
suffix = [condition,'/',subj,'/Preprocessed'];
datapath = [prefix,suffix];

%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,['*',num2str(fmod),'*preProcessed.mat'])).name};
load(datafile{1});
fname_out = [datafile{1}(1:end-4),'_matlab.mat'];
cd(cwd);
%% Data analysis & plotting:
fs = double(fs);

frames_shift = round((t_win-min(time))*fs); %shift past baseline period in MNE

pos = all_epochs_pos(:,frames_shift(1):frames_shift(2))'*1e6; %+ polarity
neg = all_epochs_neg(:,frames_shift(1):frames_shift(2))'*1e6; %- polarity

t = t_win(1):1/fs:t_win(2);
%% Get PLV spectra/Time domain waveform:

%params for random sampling with replacement
subset = 100;
k_iters = 30;

%only output things we want to look at
[f, ~, ~, PLV_env, ~, ~, T_env] = helper.getSpectAverage(pos,neg, fs, subset, k_iters);
% t = (1:length(T_env))/fs;

%% Get Peaks

[PKS,LOCS] = helper.getPeaks(f,PLV_env,fmod,harmonics);

%% Plot:
blck = [0.25, 0.25, 0.25];
rd = [0.8500, 0.3250, 0.0980, 0.5];
figure;

%Spectral Domain
hold on;
title([subj,' | RAM - 25% Duty Cycle | F_{mod}',num2str(fmod),'|',condition],'FontSize',14);
plot(f,PLV_env,'Color',blck,'linewidth',1.5);
plot(LOCS,PKS,'*','Color',rd,'MarkerSize',10,'LineWidth',2);

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
plot(t, T_env,'Color',blck, 'LineWidth',2);
xlim([0.3,.4]);
ylim([-1,1]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

set(gcf,'Position',[1557 538 560 420])

%% Export:
cd(datapath);
fname = [subj,'_RAM_',num2str(fmod),'_efr_human_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'t','T_env','f','PLV_env','PKS','LOCS')
cd(cwd);