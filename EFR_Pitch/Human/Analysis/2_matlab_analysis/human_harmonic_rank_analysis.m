clear
close all

%TODO:
% -Define colors more consistently
% -Make a more informative color scheme

%opengl a little faster/less resource heavy on linux laptop/home desktop. 
%Painters looks nicer, but either way saving figures will look fine!

set(0,'DefaultFigureRenderer','opengl')
%set(0,'DefaultFigureRenderer','painters')

tic
condition = "YNH";
type = "EFR_Pitch";
subject_no = "S353";
export = 1;
%1 to save
figsave = 1;
isPilot = 0;
local = 0;

subset = 200;
k_iters = 100;
nf_iters = 100;
F0 = 103;

%% Importing and Processing
% cwd = pwd();
% addpath(cwd);

% data_dir = strcat('../../../Data/',type,'/SNAPLab/',condition,'/',subject_no,'/FFR_Export');
% cd(data_dir);
% 
% 

cwd = pwd;
addpath(cwd);

spec = '/SNAPLab/';

if ispc && ~local
    prefix = 'A:/Pitch_Study/Data/';
elseif ~local
    [~,uname] = unix('whoami');
    uname = uname(1:end-1);
    
    if isPilot
        fold = 'Pilot_Data/';
    else
        fold = 'Pitch_Diagnostics_SH_AS/';
        spec = '/Human/';
    end
    
    prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/',fold];
else
    prefix = '/mnt/20D22780D22758F4/Shared/Code/pitch_tools/Data/';
end 

suffix = [char(type),char(spec),char(condition),'/',char(subject_no),'/EFR_Preprocessed/'];
data_dir = [prefix,suffix];

cd(data_dir);
filenames = ls(strcat('*',num2str(F0),'_alt_rank*harms_*'));

%-1 for random space at the end
filenames = split(filenames(1:(end-1))); 
ranks = zeros(length(filenames),1);
dur = .200;

for i = 1:length(filenames)
    
    file = filenames{i};
    %find rank
    ind1 = strfind(file,'rank_')+5;
    ind2 = strfind(file,'_tot_')-1;
    ranks(i) = str2double(file(ind1:ind2));
    
    load(file);
    fs = double(SampleRate); %sample rate to resample to
    onset = -window(1);
    %epoch_samps = dur*fs;
 
%     all_pos = vertcat(all_pos_fz,all_pos_cz).*1e6;
%     all_neg = vertcat(all_neg_fz,all_neg_cz)*1e6;
    
%     all_pos = 1e6*(all_pos_fz + all_pos_cz)/2;
%     all_neg = 1e6*(all_neg_fz + all_neg_cz)/2;

%       all_pos = squeeze(all_pos_fz).*1e6;
%       all_neg = squeeze(all_neg_fz).*1e6;
% 
    all_pos = squeeze(mean(all_pos_cz,2))*1e6;
    all_neg = squeeze(mean(all_neg_cz,2))*1e6;
% 
%     all_pos = squeeze(mean(all_pos_central,2))*1e6;
%     all_neg = squeeze(mean(all_neg_central,2))*1e6;

    pos_abr = all_pos';
    neg_abr = all_neg';

    trials = size(pos_abr,2);
    nsamps = size(pos_abr,1);
            
    s1 = round(onset*fs);
    s2 = round((onset+dur)*fs);
    
    [l_floorx,l_floory] = helper.getNoiseFloor(pos_abr(s1:s2,:),neg_abr(s1:s2,:),fs,nf_iters);
    
    disp('Getting Spectral Average')
    
    [f, MAG_env, MAG_tfs, PLV_env, PLV_tfs, T_tfs, T_env] = helper.getSpectAverage(pos_abr(s1:s2,:),neg_abr(s1:s2,:), fs, subset, k_iters);
    
    %Hacky....change this later
    PLV_env_a_mean(:,i) = PLV_env;
    PLV_tfs_a_mean(:,i) = PLV_tfs;
    
    T_env_a_mean(:,i) = T_env;
    T_tfs_a_mean(:,i) = T_tfs;
    
    MAG_env_a_mean(:,i) = MAG_env;
    a_floor_mean(:,i) = l_floory;
end

mag_env = db2mag(abs(MAG_env_a_mean))-db2mag(a_floor_mean);
phase_env = angle(MAG_env_a_mean);

est_harmonics = 4;
for p = 1:length(filenames)
    [amp_harms_env, loc] = helper.getPeaks(f,mag_env(:,p),F0,est_harmonics); 
    t = 0:1/fs:dur;
    for h = 1:est_harmonics
        est_wform(:,h) = amp_harms_env(h).*sin(2*pi*h*F0.*t+phase_env(loc(h),p)); 
    end
    est_combined = sum(est_wform,2);
    T_env_est(:,p) = est_combined/max(abs(est_combined));
end

T_env_est_mean = T_env_est;

% cd(cwd)

%% Sorting matrices by ascending rank
[ranks_ord,I] = sort(ranks,'ascend');
T_env_a_mean = T_env_a_mean(:,I);
T_env_est_mean = T_env_est_mean(:,I);
PLV_env_a_mean = PLV_env_a_mean(:,I);
PLV_tfs_a_mean = PLV_tfs_a_mean(:,I);

%% Make EFR Figure:
efrfig = figure;
subplot(1,3,1);
buff = 0.3;
time = (1:length(T_tfs_a_mean))/fs;
%Filter EFRs
[b,a] = butter(6,[20,2000]/(fs/2));
T_env_unfilt = T_env_a_mean;
T_env_a_mean = filtfilt(b,a,T_env_a_mean);

hold on;
for r = 1:length(ranks_ord)
    plot(time,T_env_a_mean(:,r)+buff*(ranks_ord(r)-1),'linewidth',1.5)
end

xlim([0,max(time)])
hold off
xlabel('Time (s)')
title('Envelope-Following Response')
ylabel('Amplitude (\muV)')
yticks([]);
rectangle('Position',[0.01,0.15,0.001,.5],'facecolor','k')
text(0.012,0.3,['.5 \muV'],'FontWeight','Bold');
%% Process using windowed autocorelation:

window = 1.5/F0; %window in seconds
acf_bin_length = round(window*fs);
inds = 1:acf_bin_length/10:(length(T_tfs_a_mean)-round(window*fs));
inds = round(inds);

%Running ACF
[acf_mean,t_acf,P1,P2] = helper.run_ACF_Pitch(T_env_a_mean,F0,fs,0);

figure(efrfig);
subplot(1,3,2);
plot(t_acf,acf_mean,'LineWidth',2);
xlabel('Time-lag(ms)');

text(1000/(2*F0),.6,'P2')
text(1000/F0,.4,'P1')
xlabel('Time-lag (ms)')
legend(num2str(ranks_ord),'location','NorthEast');
[hleg,att] = legend('show');
title(hleg,'Harmonic Rank')
title('Auto-Correlation Function')
ylabel('Normalized Amplitude')

%% Sigmoid

%calculate strength of 2nd Harmonic
sig_data = P2./(P1+P2);
% sig_data = P2./P1;


%Fit model sigmoid
x = 0:0.1:15;
maximum = 1.2;
mid =6;
steep = 1.3;
start = 0.01;
sigmoid = 'a./(1+exp(-b*(x-c)))+d';
startPoints = [maximum, steep, mid, start];
fops = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf, 0, 1, 0],'Upper',[inf, inf, 15, inf],'StartPoint',startPoints);
ft = fittype(sigmoid,'options',fops);

sig_fit = fit(ranks_ord, sig_data,ft);
sig_model = sig_fit(x);

subplot(1,3,3)
hold on
for k = 1:length(ranks_ord)
    scatter(ranks_ord(k),sig_data(k),50,'filled');
end 
plot(x,sig_model,'--k','LineWidth',2);
xlabel('Harmonic Rank')
title('Strength of 2nd Harmonic ACF Ratio')
ylabel('(P2/(P1+P2)')
ylim([0,1])
xlim([1,12])
xticks(ranks_ord)
xticklabels(ranks_ord)
hold off
efrtitle = strcat(subject_no,' | ',condition);
sgtitle(efrtitle)
%% Spectral Envelope Analysis:

plvfig = figure;
subplot(1,2,1)
hold on
buff = 0.07;

for k = 1:length(ranks_ord)
    plot(f,PLV_env_a_mean(:,k)+buff*(ranks_ord(k)),'linewidth',1.5);
end

title('ENV');
xlabel('Frequency (Hz)')
ylabel('PLV');
yticks([]);
xlim([0,4000])
ylim([buff,(max(ranks_ord)+1.5)*buff]);
%% Spectral TFS Analysis:

subplot(1,2,2);
hold on
buff = 0.07;

for k = 1:length(ranks_ord)
    plot(f,PLV_tfs_a_mean(:,k)+buff*(ranks_ord(k)),'linewidth',1.5);
end

title('TFS');
xlabel('Frequency (Hz)')
legend(num2str(ranks_ord),'location','NorthEast');
[hleg,att] = legend('show');
title(hleg,'Harmonic Rank')
yticks([]);
ylim([buff,(max(ranks_ord)+1.5)*buff]);
yticklabels('');
xlim([0,4000])
plvtitle = strcat(subject_no,' | ',condition,' - PLV Spectra');
sgtitle(plvtitle)

% %% Estimated EFR figure
% estefr = figure();
% subplot(1,3,1);
% [ranks_ord,I] = sort(ranks,'ascend');
% buff = 1;
% time = (1:length(T_tfs_a_mean))/fs;
% %Filter EFRs
% [b,a] = butter(6,[65,400]/(fs/2));
% T_env_a_mean = filtfilt(b,a,T_env_a_mean);
% 
% hold on;
% for r = 1:length(ranks_ord)
%     plot(time,T_env_est(:,r)+buff*(ranks_ord(r)-1),'linewidth',1.5)
% 
% end
% xlim([0,max(time)])
% hold off
% xlabel('Time (s)')
% title('Estimated Envelope-Following Response')
% ylabel('Amplitude (Arbitrary)')
% 
% % ACF
% [est_acf_mean,t_acf,est_P1,est_P2] = run_ACF_Pitch(T_env_est_mean,F0,fs);
% 
% figure(estefr)
% subplot(1,3,2);
% plot(t_acf,est_acf_mean,'LineWidth',2);
% xlabel('Time-lag(ms)');
% text(1000/(2*F0),.6,'P2')
% text(1000/F0,.4,'P1')
% xlabel('Time-lag (ms)')
% legend(num2str(ranks_ord),'location','NorthEast');
% [hleg,att] = legend('show');
% title(hleg,'Harmonic Rank')
% title('Auto-Correlation Function')
% ylabel('Normalized Amplitude')
% 
% %Est Sigmoid
% est_sig_data = est_P2./(est_P1+est_P2);
% %Fit model sigmoid
% x = 0:0.1:15;
% maximum = 1.2;
% mid =6;
% steep = 1.3;
% start = 0.01;
% sigmoid = 'a./(1+exp(-b*(x-c)))+d';
% startPoints = [maximum, steep, mid, start];
% fops = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf, 0, 1, 0],'Upper',[inf, inf, 15, inf],'StartPoint',startPoints);
% ft = fittype(sigmoid,'options',fops);
% 
% sig_fit = fit(ranks_ord, est_sig_data,ft);
% sig_model = sig_fit(x);
% 
% subplot(1,3,3)
% hold on
% for k = 1:length(ranks_ord)
%     scatter(ranks_ord(k),sig_data(k),50,'filled');
% end 
% plot(x,sig_model,'--k','LineWidth',2);
% xlabel('Harmonic Rank')
% title('Strength of 2nd Harmonic ACF Ratio')
% ylabel('(P2/(P1+P2)')
% ylim([0,1])
% xlim([3,17])
% xticks(ranks_ord)
% xticklabels(ranks_ord)
% hold off
% efrtitle = strcat(subject_no,' | Estimated EFR ',condition);
% sgtitle(efrtitle)
%% Save

    set(efrfig, 'Position', [192,401,1510,553]);
%     set(estefr, 'Position', [192,401,1510,553]);
    set(plvfig, 'Position', [887,26,1014,501]);

if(figsave~=0)
%     if(~exist('Processed'))
%         mkdir('Processed');
%     end
% 
%     cd('Processed')
    cd(data_dir);
    print(efrfig,strcat(subject_no,'_',condition,'_EFR_processed'),'-r600','-dpng')
    print(plvfig,strcat(subject_no,'_',condition,'_PLV_processed'),'-r600','-dpng')
%     print(acf_fig,strcat(subject_no,'_',condition,'_ACG_processed'),'-r600','-dpng')
%     print(estefr,strcat(subject_no,'_',condition,'_EFR_estimated'),'-r600','-dpng')
    cd ../
end

%% export mat file
if(export~=0)
%     if(~exist('Processed'))
%         mkdir('Processed');
%     end
%     cd('Processed')
    cd(data_dir);
    save(strcat(subject_no,'_',condition,'_data.mat'),'time', 'T_env_unfilt','T_env_a_mean','t_acf','acf_mean','P1','P2','ranks_ord','PLV_tfs_a_mean','PLV_env_a_mean','f');
end

cd(cwd)

%% Only Time-Domain Waveform:

% wfmfig = figure;
% hold on;
% 
% buff = 0.3;
% factor = 0.1;
% c_plot = [0 0.3 0];
% for r = 1:length(ranks_ord)
%     plot(time,T_env_a_mean(:,r)+buff*(ranks_ord(r)-1),'linewidth',1.5,'Color',c_plot.^factor)
%     text(0.06,buff*(ranks_ord(r)-1)-0.4,strcat("Rank = ",num2str(ranks_ord(r))),'Color',c_plot.^factor,'fontweight','bold');
%     factor = factor*(1+i*factor);
% end
% 
% xlim([0,.200])
% hold off
% xlabel('Time (s)')
% title('Envelope-Following Response')
% ylabel('Amplitude (\muV)')
% yticks([]);
% rectangle('Position',[0.005,0.15,0.001,.5],'facecolor','k')
% text(0.007,0.3,['.5 \muV'],'FontWeight','Bold');
% 
% cd('Processed')
% print(wfmfig,strcat(subject_no,'_',condition,'human_efr'),'-dsvg')
% cd(cwd)
