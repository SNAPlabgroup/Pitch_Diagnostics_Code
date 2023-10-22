%Author: Andrew Sivaprakasam
%Updated: July, 2023
%Purpose: Script to import/plot/apply additional processing to RAM_EFR
%files
%Helpful Info: Be sure to define datapath so Import data section works. 
%see my example setup_AS file.   

%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,[subj,'_ABR_preProcessed.mat'])).name};
load(datafile{1});
fname_out = [datafile{1}(1:end-4),'_matlab.mat'];
cd(cwd);
%% Peak Finding
% local max / peaks
peak = islocalmax(abr);
valley = islocalmin(abr); 
peak_new = zeros(size(peak)); 
valley_new = zeros(size(valley)); 
for n = 1:length(t)
    if (t(n) > 1.9) && (t(n) < 6.2)
        if peak(n)
            peak_new(n) = 1; 
        end
    end
end
for n = 1:length(t)
    if (t(n) > 2.4) && (t(n) < 7.2)
        if valley(n)
            valley_new(n) = 1; 
        end
    end
end
peak_new=logical(peak_new); 
valley_new=logical(valley_new); 

peak_time = t(peak_new); 
valley_time = t(valley_new); 
peak_amp = abr(peak_new); 
valley_amp = abr(valley_new); 

%% Wave Amp
w1_amp = peak_amp(1) - valley_amp(1); 
w5_amp = peak_amp(5) - valley_amp(5);
I_V_ratio = w1_amp/w5_amp; 
%% Plotting
fs = double(fs);
xstart = .6;
xend = .9;
ystart = 0.6;
yend = .9;
figure; 
t = time_s*1e3; % to ms
trials = min(size(all_epochs_neg, 1),size(all_epochs_pos, 1)); 
abr = combined_mean.*1e6; % to uV; 

%axes('Position',[xstart ystart xend-xstart yend-ystart])
box on
hold on
plot(t, abr,'Color','k', 'LineWidth',2);
xlim([0,10]);
ylim([-0.5,0.5]);
yticks([-.5:.1:.5])
xlabel('Time (ms)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
title('Auditory Brainstem Response')

text(peak_time(1)-.035, peak_amp(1)+0.03, 'I', 'FontSize', 14)
text(peak_time(2)-0.06, peak_amp(2)+0.03, 'II', 'FontSize', 14)
text(peak_time(3)-0.1, peak_amp(3)+0.03, 'III', 'FontSize', 14)
text(peak_time(4)-0.1, peak_amp(4)+0.03, 'IV', 'FontSize', 14)
text(peak_time(5)-0.09, peak_amp(5)+0.03, 'V', 'FontSize', 14)

%plot(peak_time, peak_amp+.03,'rv', 'MarkerSize', 10, 'linew', 1.5)
plot([valley_time(1), valley_time(5)], [valley_amp(1)-.03, valley_amp(5)-.03],'b^', 'MarkerSize', 10, 'linew', 1.5)
%plot([peak_time(1), peak_time(5)], [peak_amp(1)-.03, peak_amp(5)-.03],'r^', 'MarkerSize', 10, 'linew', 1.5)

set(gcf,'Position',[1557 538 560 420]); 
set(gca,'FontName', 'Times', 'FontSize', 16)

%% Export:
cd(datapath);
fname = [subj,'_abr_human_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'valley_amp', 'valley_time', 'peak_amp', 'peak_time', 'I_V_ratio','abr','t' )
cd(cwd);


% %% Data analysis & plotting:
% 
% 
% frames_shift = round((t_win-min(time))*fs); %shift past baseline period in MNE
% 
% pos = all_epochs_pos(:,frames_shift(1):frames_shift(2))'*1e6; %+ polarity
% neg = all_epochs_neg(:,frames_shift(1):frames_shift(2))'*1e6; %- polarity
% 
% t = t_win(1):1/fs:t_win(2);
% %% Get PLV spectra/Time domain waveform:
% 
% %params for random sampling with replacement
% subset = 100;
% k_iters = 30;
% 
% %only output things we want to look at
% [f, ~, ~, PLV_env, ~, ~, T_env] = helper.getSpectAverage(pos,neg, fs, subset, k_iters);
% % t = (1:length(T_env))/fs;
% 
% %% Get Peaks
% 
% [PKS,LOCS] = helper.getPeaks(f,PLV_env,fmod,harmonics);
% 
% %% Plot:
% blck = [0.25, 0.25, 0.25];
% rd = [0.8500, 0.3250, 0.0980, 0.5];
% figure;
% 
% %Spectral Domain
% hold on;
% title([subj,' | RAM - 25% Duty Cycle | F_{mod}',num2str(fmod),'|',condition],'FontSize',14);
% plot(f,PLV_env,'Color',blck,'linewidth',1.5);
% plot(LOCS,PKS,'*','Color',rd,'MarkerSize',10,'LineWidth',2);
% 
% hold off;
% ylim([0,1])
% ylabel('PLV','FontWeight','bold')
% xlabel('Frequency(Hz)','FontWeight','bold')
% 
% %Time Domain