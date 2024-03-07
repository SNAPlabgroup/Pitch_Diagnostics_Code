%% Script to plot all data collected

clear;
chins = {'Q412', 'Q423', 'Q424', 'Q426', 'Q430', 'Q431', 'Q427', 'Q428', 'Q421', 'Q425', 'Q443', 'Q422', 'Q440', 'Q441'};
group = {'TTS', 'TTS', 'TTS', 'TTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA'};

location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

TTS_count = 0;
PTS_count = 0;
CA_count = 0;
GE_count = 0;

codedir = pwd; 

awake = 1; 

%% GET DATA
for k = 1:length(chins)
    
    subj = chins{k};
    
    if strcmp(group(k), 'TTS')
        TTS_count = TTS_count+1; %add new chin to the list
        TTS_chin{TTS_count} = chins{k};
        
        
        % first get baseline data
        [TTS_Pre_oae(TTS_count, :),TTS_Pre_nf(TTS_count,:), ...
            TTS_Pre_f2(TTS_count,:),TTS_Pre_oaesum(TTS_count,:), ...
            TTS_Pre_centerfreq(TTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 1day Post data
        [TTS_Post_1d_oae(TTS_count, :),TTS_Post_1d_nf(TTS_count,:), ...
            TTS_Post_1d_f2(TTS_count,:),TTS_Post_1d_oaesum(TTS_count,:), ...
            TTS_Post_1d_centerfreq(TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_1dayPost', awake);
        
        % Then get 2wks Post data
        [TTS_Post_oae(TTS_count, :),TTS_Post_nf(TTS_count,:), ...
            TTS_Post_f2(TTS_count,:),TTS_Post_oaesum(TTS_count,:), ...
            TTS_Post_centerfreq(TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_2wksPost', awake);
        
    elseif strcmp(group{k}, 'PTS')
        PTS_count = PTS_count+1; %add new chin to the list
        PTS_chin{PTS_count} = chins{k};
        
        % first get baseline data
        [PTS_Pre_oae(PTS_count, :),PTS_Pre_nf(PTS_count,:), ...
            PTS_Pre_f2(PTS_count,:),PTS_Pre_oaesum(PTS_count,:), ...
            PTS_Pre_centerfreq(PTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 2wks Post data
        [PTS_Post_oae(PTS_count, :),PTS_Post_nf(PTS_count,:), ...
            PTS_Post_f2(PTS_count,:),PTS_Post_oaesum(PTS_count,:), ...
            PTS_Post_centerfreq(PTS_count,:)] = get_data_DP(prefix, subj, 'PTS_2wksPost', awake);
        
        
    elseif strcmp(group{k}, 'CA')
        CA_count = CA_count+1; %add new chin to the list
        CA_chin{CA_count} = chins{k};
        
        % first get baseline data
        [CA_Pre_oae(CA_count, :),CA_Pre_nf(CA_count,:), ...
            CA_Pre_f2(CA_count,:),CA_Pre_oaesum(CA_count,:), ...
            CA_Pre_centerfreq(CA_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 2wks Post data
        [CA_Post_oae(CA_count, :),CA_Post_nf(CA_count,:), ...
            CA_Post_f2(CA_count,:),CA_Post_oaesum(CA_count,:), ...
            CA_Post_centerfreq(CA_count,:)] = get_data_DP(prefix, subj, 'CA_2wksPost', awake);
        
        
    elseif strcmp(group{k}, 'GE')
        GE_count = GE_count+1; %add new chin to the list
        GE_chin{GE_count} = chins{k};
        
        % first get baseline data
        [GE_Pre_oae(GE_count, :),GE_Pre_nf(GE_count,:), ...
            GE_Pre_f2(GE_count,:),GE_Pre_oaesum(GE_count,:), ...
            GE_Pre_centerfreq(GE_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 2wks Post data
        [GE_Post_oae(GE_count, :),GE_Post_nf(GE_count,:), ...
            GE_Post_f2(GE_count,:),GE_Post_oaesum(GE_count,:), ...
            GE_Post_centerfreq(GE_count,:)] = get_data_DP(prefix, subj, 'GE_2wksPost', awake);
        
        
    end
    
end

%% Get averages

%% Plot Data

blck = [0.25, 0.25, 0.25];
rd = [217, 95, 2]./255; %TTS
blu = [117, 112, 179]./255; %CA
yel = [27, 158, 119]./255; %PTS
gre = [115, 177, 117]./255; %GE

i_blck = [0.75, 0.75, 0.75];
i_rd = [217, 95, 2 75]./255; %TTS
i_blu = [117, 112, 179, 75]./255; %CA
i_yel = [27, 158, 119, 75]./255; %PTS
i_gre = [115, 177, 117, 57]./255; %GE

figure; 
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])
% TTS
subplot(2,2,1)
hold on;
%text(.71,45, sprintf('n = %d', size(TTS_Pre_f2, 1) - sum(isnan(TTS_Pre_f2(:,1)))), 'FontSize', 12)
%text(.71,40, sprintf('n = %d', size(TTS_Post_f2, 1) - sum(isnan(TTS_Post_f2(:,1)))), 'Color', rd, 'FontSize', 12)
%text(.71,35, sprintf('n = %d', size(TTS_Post_1d_f2, 1) - sum(isnan(TTS_Post_1d_f2(:,1)))), 'Color', rd, 'FontSize', 12)
% for i = 1:TTS_count
%     plot(TTS_Pre_f2(i,:), TTS_Pre_oae(i,:),'Color',i_blck,'linewidth',3);
%     plot(TTS_Post_1d_f2(i,:), TTS_Post_1d_oae(i,:),':', 'Color', i_rd, 'linew', 3); 
%     plot(TTS_Post_f2(i,:), TTS_Post_oae(i,:), 'Color', i_rd, 'linew', 3); 
% end

plot(mean(TTS_Pre_f2,1, 'omitNaN'),mean(TTS_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',4);
plot(mean(TTS_Pre_f2,1, 'omitNaN'),mean(TTS_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',4);
plot(mean(TTS_Pre_centerfreq,1, 'omitNaN'),mean(TTS_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',4);

hold on;
plot(mean(TTS_Post_1d_f2,1, 'omitNaN'),mean(TTS_Post_1d_oae,1, 'omitNaN'),':', 'Color',rd,'linewidth',4);
plot(mean(TTS_Post_1d_f2,1, 'omitNaN'),mean(TTS_Post_1d_nf,1, 'omitNaN'),':', 'Color',rd,'linewidth',4);
plot(mean(TTS_Post_1d_centerfreq,1, 'omitNaN'),mean(TTS_Post_1d_oaesum,1, 'omitNaN'),'*','Color',rd,'MarkerSize',10,'LineWidth',4);

hold on;
plot(mean(TTS_Post_f2,1, 'omitNaN'),mean(TTS_Post_oae,1, 'omitNaN'),'Color',rd,'linewidth',4);
plot(mean(TTS_Post_f2,1, 'omitNaN'),mean(TTS_Post_nf,1, 'omitNaN'),'--', 'Color',rd,'linewidth',4);
plot(mean(TTS_Post_centerfreq,1, 'omitNaN'),mean(TTS_Post_oaesum,1, 'omitNaN'),'*','Color',rd,'MarkerSize',10,'LineWidth',4);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency (kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 20)
title('Synaptopathy','FontSize',24, 'Color', rd);
grid on; 

% PTS
subplot(2,2,2)
hold on;
%text(.71,45, sprintf('n = %d', size(PTS_Pre_f2, 1) - sum(isnan(PTS_Pre_f2(:,1)))), 'FontSize', 12)
%text(.71,40, sprintf('n = %d', size(PTS_Post_f2, 1) - sum(isnan(PTS_Post_f2(:,1)))), 'Color', yel, 'FontSize', 12)
% for i = 1:PTS_count
%     plot(PTS_Pre_f2(i,:), PTS_Pre_oae(i,:),'Color',i_blck,'linewidth',3);
%     plot(PTS_Post_f2(i,:), PTS_Post_oae(i,:), 'Color', i_yel, 'linew', 3); 
% end

plot(mean(PTS_Pre_f2,1, 'omitNaN'),mean(PTS_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',4);
plot(mean(PTS_Pre_f2,1, 'omitNaN'),mean(PTS_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',4);
plot(mean(PTS_Pre_centerfreq,1, 'omitNaN'),mean(PTS_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',4);

hold on;
plot(mean(PTS_Post_f2,1, 'omitNaN'),mean(PTS_Post_oae,1, 'omitNaN'),'Color',yel,'linewidth',4);
plot(mean(PTS_Post_f2,1, 'omitNaN'),mean(PTS_Post_nf,1, 'omitNaN'),'--', 'Color',yel,'linewidth',4);
plot(mean(PTS_Post_centerfreq,1, 'omitNaN'),mean(PTS_Post_oaesum,1, 'omitNaN'),'*','Color',yel,'MarkerSize',10,'LineWidth',4);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency (kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 20)
title('Complex SNHL','FontSize',24, 'Color', yel);
grid on; 

% CA
subplot(2,2,3)
hold on;
%text(.71,45, sprintf('n = %d', size(CA_Pre_f2, 1) - sum(isnan(CA_Pre_f2(:,1)))), 'FontSize', 12)
%text(.71,40, sprintf('n = %d', size(CA_Post_f2, 1) - sum(isnan(CA_Post_f2(:,1)))), 'Color', blu, 'FontSize', 12)
hold on; 
% for i = 1:CA_count
%     plot(CA_Pre_f2(i,:), CA_Pre_oae(i,:),'Color',i_blck,'linewidth',3);
%     plot(CA_Post_f2(i,:), CA_Post_oae(i,:), 'Color', i_blu, 'linew', 3); 
% end
plot(mean(CA_Pre_f2,1, 'omitNaN'),mean(CA_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',4);
plot(mean(CA_Pre_f2,1, 'omitNaN'),mean(CA_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',4);
plot(mean(CA_Pre_centerfreq,1, 'omitNaN'),mean(CA_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',4);

hold on;
plot(mean(CA_Post_f2,1, 'omitNaN'),mean(CA_Post_oae,1, 'omitNaN'),'Color',blu,'linewidth',4);
plot(mean(CA_Post_f2,1, 'omitNaN'),mean(CA_Post_nf,1, 'omitNaN'),'--', 'Color',blu,'linewidth',4);
plot(mean(CA_Post_centerfreq,1, 'omitNaN'),mean(CA_Post_oaesum,1, 'omitNaN'),'*','Color',blu,'MarkerSize',10,'LineWidth',4);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency (kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 20)
title('IHC Dysfunction', 'FontSize', 24, 'Color', blu)
grid on; 

% subplot(2,2,4)
% hold on;
% title('GE | Swept DPOAEs','FontSize',16);
% text(.71,45, sprintf('n = %d', size(GE_Pre_f2, 1) - sum(isnan(GE_Pre_f2(:,1)))), 'FontSize', 12)
% text(.71,40, sprintf('n = %d', size(GE_Post_f2, 1) - sum(isnan(GE_Post_f2(:,1)))), 'Color', gre, 'FontSize', 12)
% plot(mean(GE_Pre_f2,1, 'omitNaN'),mean(GE_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',3);
% plot(mean(GE_Pre_f2,1, 'omitNaN'),mean(GE_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',3);
% plot(mean(GE_Pre_centerfreq,1, 'omitNaN'),mean(GE_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',3);
% 
% hold on;
% plot(mean(GE_Post_f2,1, 'omitNaN'),mean(GE_Post_oae,1, 'omitNaN'),'Color',gre,'linewidth',3);
% plot(mean(GE_Post_f2,1, 'omitNaN'),mean(GE_Post_nf,1, 'omitNaN'),'--', 'Color',gre,'linewidth',3);
% plot(mean(GE_Post_centerfreq,1, 'omitNaN'),mean(GE_Post_oaesum,1, 'omitNaN'),'*','Color',gre,'MarkerSize',10,'LineWidth',3);
% hold off;
% 
% ylim([-40, 50])
% xlim([.5, 16])
% xticks([.5, 1, 2, 4, 8, 16])
% ylabel('Amplitude (dB EPL)','FontWeight','bold')
% xlabel('Frequency(kHz)','FontWeight','bold')
% set(gca, 'XScale', 'log', 'FontSize', 14)

cd 'Figures'
print -dpng -r600 Pre-Post-DP
cd .. 
