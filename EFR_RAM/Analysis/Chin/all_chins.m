%% Script to plot all data collected

clear;

chins = {'Q412', 'Q424', 'Q426', 'Q430', 'Q431', 'Q440', 'Q441', 'Q427', 'Q428', 'Q425', 'Q421', 'Q422'};
group = {'TTS', 'TTS', 'TTS', 'CA', 'CA', 'CA', 'CA','PTS', 'PTS', 'CA', 'CA',  'PTS'};

conditions = {'Baseline', '_2wksPost'};

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

%% GET DATA
for k = 1:length(chins)
    
    subj = chins{k};
    
    if strcmp(group(k), 'TTS')
        TTS_count = TTS_count+1; %add new chin to the list
        TTS_chin{TTS_count} = chins{k};
        
        % first get baseline data
        [TTS_Pre_t_all(TTS_count, :),TTS_Pre_T_env_all(TTS_count,:), ...
            TTS_Pre_PKS_all(TTS_count,:),TTS_Pre_PLV_env_all(TTS_count,:), ...
            TTS_Pre_f_all(TTS_count,:), TTS_Pre_LOCS_all(TTS_count,:)] = get_data(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [TTS_Post_t_all(TTS_count, :),TTS_Post_T_env_all(TTS_count,:), ...
            TTS_Post_PKS_all(TTS_count,:),TTS_Post_PLV_env_all(TTS_count,:), ...
            TTS_Post_f_all(TTS_count,:), TTS_Post_LOCS_all(TTS_count,:)] = get_data(prefix, subj, 'TTS_2wksPost');
        
    elseif strcmp(group{k}, 'PTS')
        PTS_count = PTS_count+1; %add new chin to the list
        PTS_chin{PTS_count} = chins{k};
        
        % first get baseline data
        [PTS_Pre_t_all(PTS_count, :),PTS_Pre_T_env_all(PTS_count,:), ...
            PTS_Pre_PKS_all(PTS_count,:),PTS_Pre_PLV_env_all(PTS_count,:), ...
            PTS_Pre_f_all(PTS_count,:), PTS_Pre_LOCS_all(PTS_count,:)] = get_data(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [PTS_Post_t_all(PTS_count, :),PTS_Post_T_env_all(PTS_count,:), ...
            PTS_Post_PKS_all(PTS_count,:),PTS_Post_PLV_env_all(PTS_count,:), ...
            PTS_Post_f_all(PTS_count,:), PTS_Post_LOCS_all(PTS_count,:)] = get_data(prefix, subj, 'PTS_2wksPost');
        
        
    elseif strcmp(group{k}, 'CA')
        CA_count = CA_count+1; %add new chin to the list
        CA_chin{CA_count} = chins{k};
        
        % first get baseline data
        [CA_Pre_t_all(CA_count, :),CA_Pre_T_env_all(CA_count,:), ...
            CA_Pre_PKS_all(CA_count,:),CA_Pre_PLV_env_all(CA_count,:), ...
            CA_Pre_f_all(CA_count,:), CA_Pre_LOCS_all(CA_count,:)] = get_data(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [CA_Post_t_all(CA_count, :),CA_Post_T_env_all(CA_count,:), ...
            CA_Post_PKS_all(CA_count,:),CA_Post_PLV_env_all(CA_count,:), ...
            CA_Post_f_all(CA_count,:), CA_Post_LOCS_all(CA_count,:)] = get_data(prefix, subj, 'CA_2wksPost');
        
        
    elseif strcmp(group{k}, 'GE')
        GE_count = GE_count+1; %add new chin to the list
        GE_chin{GE_count} = chins{k};
        
        % first get baseline data
        [GE_Pre_t_all(GE_count, :),GE_Pre_T_env_all(GE_count,:), ...
            GE_Pre_PKS_all(GE_count,:),GE_Pre_PLV_env_all(GE_count,:), ...
            GE_Pre_f_all(GE_count,:), GE_Pre_LOCS_all(GE_count,:)] = get_data(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [GE_Post_t_all(GE_count, :),GE_Post_T_env_all(GE_count,:), ...
            GE_Post_PKS_all(GE_count,:),GE_Post_PLV_env_all(GE_count,:), ...
            GE_Post_f_all(GE_count,:), GE_Post_LOCS_all(GE_count,:)] = get_data(prefix, subj, 'GE_2wksPost');
        
        
    end
    
end

%% Get ratio
% CA_pre_ratio = sum(CA_Pre_PKS_all(:, 3:8),2)./sum(CA_Pre_PKS_all(:,1:2),2);
% CA_post_ratio = sum(CA_Post_PKS_all(:, 3:8),2)./sum(CA_Post_PKS_all(:,1:2),2);
% PTS_pre_ratio = sum(PTS_Pre_PKS_all(:, 3:8),2)./sum(PTS_Pre_PKS_all(:,1:2),2);
% PTS_post_ratio = sum(PTS_Post_PKS_all(:, 3:8),2)./sum(PTS_Post_PKS_all(:,1:2),2);
% TTS_pre_ratio = sum(TTS_Pre_PKS_all(:, 3:8),2)./sum(TTS_Pre_PKS_all(:,1:2),2);
% TTS_post_ratio = sum(TTS_Post_PKS_all(:, 3:8),2)./sum(TTS_Post_PKS_all(:,1:2),2);


harms_locs(1) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,1)); 
harms_locs(2) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,2)); 
harms_locs(3) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,3)); 
harms_locs(4) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,4)); 
harms_locs(5) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,5)); 
% harms_locs(6) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,6)); 
% harms_locs(7) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,7)); 
% harms_locs(8) = find(PTS_Pre_f_all(1,:) == PTS_Pre_LOCS_all(1,8)); 



CA_pre_ratio = sum(CA_Pre_PLV_env_all(:, harms_locs),2); %./sum(CA_Pre_PKS_all(:,1:3),2);
CA_post_ratio = sum(CA_Post_PLV_env_all(:, harms_locs),2); %./sum(CA_Post_PKS_all(:,1:3),2);
PTS_pre_ratio = sum(PTS_Pre_PLV_env_all(:, harms_locs),2); %./sum(PTS_Pre_PKS_all(:,1:3),2);
PTS_post_ratio = sum(PTS_Post_PLV_env_all(:, harms_locs),2); %./sum(PTS_Post_PKS_all(:,1:3),2);
TTS_pre_ratio = sum(TTS_Pre_PLV_env_all(:, harms_locs),2); %./sum(TTS_Pre_PKS_all(:,1:3),2);
TTS_post_ratio = sum(TTS_Post_PLV_env_all(:, harms_locs),2); %./sum(TTS_Post_PKS_all(:,1:3),2);

%% Plot Data

blck = [0.25, 0.25, 0.25];
rd = [217, 95, 2]./255; %TTS
blu = [117, 112, 179]./255; %CA
yel = [27, 158, 119]./255; %PTS
gre = [115, 177, 117]./255; %GE

%Time Domain
xstart = .2;
xend = .5;
ystart = 0.6;
yend = .9;

figure; %Spectral Domain
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])
% TTS
subplot(2,2,1)
hold on;
set(gca, 'FontSize', 20)
title('Synaptopathy','FontSize',24, 'Color', rd);
plot(mean(TTS_Pre_f_all,1, 'omitNaN'),mean(TTS_Pre_PLV_env_all,1, 'omitNaN'),'Color',blck,'linewidth',3);
plot(mean(TTS_Pre_LOCS_all,1, 'omitNaN'),mean(TTS_Pre_PKS_all,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',3);

hold on;
plot(mean(TTS_Post_f_all,1, 'omitNaN'),mean(TTS_Post_PLV_env_all,1, 'omitNaN'),'Color',rd,'linewidth',3);
plot(mean(TTS_Post_LOCS_all,1, 'omitNaN'),mean(TTS_Post_PKS_all,1, 'omitNaN'),'*','Color',rd,'MarkerSize',10,'LineWidth',3);

hold off;
ylim([0,1])
ylabel('PLV','FontWeight','bold')
xlabel('Frequency (Hz)','FontWeight','bold')

% Carbo
subplot(2,2,3)
hold on;
set(gca, 'FontSize', 20)
title('IHC Dysfunction','FontSize',24, 'Color', blu);
plot(mean(CA_Pre_f_all,1, 'omitNaN'),mean(CA_Pre_PLV_env_all,1, 'omitNaN'),'Color',blck,'linewidth',3);
plot(mean(CA_Pre_LOCS_all,1, 'omitNaN'),mean(CA_Pre_PKS_all,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',3);

hold on;
plot(mean(CA_Post_f_all,1, 'omitNaN'),mean(CA_Post_PLV_env_all,1, 'omitNaN'),'Color',blu,'linewidth',3);
plot(mean(CA_Post_LOCS_all,1, 'omitNaN'),mean(CA_Post_PKS_all,1, 'omitNaN'),'*','Color',blu,'MarkerSize',10,'LineWidth',3);

hold off;
ylim([0,1])
ylabel('PLV','FontWeight','bold')
xlabel('Frequency (Hz)','FontWeight','bold')

% PTS
subplot(2,2,2)
hold on;
set(gca, 'FontSize', 20)
title('Complex SNHL','FontSize',24, 'Color', yel);
plot(mean(PTS_Pre_f_all,1, 'omitNaN'),mean(PTS_Pre_PLV_env_all,1, 'omitNaN'),'Color',blck,'linewidth',3);
plot(mean(PTS_Pre_LOCS_all,1, 'omitNaN'),mean(PTS_Pre_PKS_all,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',3);

hold on;
plot(mean(PTS_Post_f_all,1, 'omitNaN'),mean(PTS_Post_PLV_env_all,1, 'omitNaN'),'Color',yel,'linewidth',3);
plot(mean(PTS_Post_LOCS_all,1, 'omitNaN'),mean(PTS_Post_PKS_all,1, 'omitNaN'),'*','Color',yel,'MarkerSize',10,'LineWidth',3);

hold off;
ylim([0,1])
ylabel('PLV','FontWeight','bold')
xlabel('Frequency (Hz)','FontWeight','bold')

% subplot(2,2,4)
% hold on;
% title('GE | RAM - 25% Duty Cycle','FontSize',14);
% plot(mean(GE_Pre_f_all,1, 'omitNaN'),mean(GE_Pre_PLV_env_all,1, 'omitNaN'),'Color',blck,'linewidth',1.5);
% plot(mean(GE_Pre_LOCS_all,1, 'omitNaN'),mean(GE_Pre_PKS_all,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',2);
% 
% hold on;
% plot(mean(GE_Post_f_all,1, 'omitNaN'),mean(GE_Post_PLV_env_all,1, 'omitNaN'),'Color',gre,'linewidth',1.5);
% plot(mean(GE_Post_LOCS_all,1, 'omitNaN'),mean(GE_Post_PKS_all,1, 'omitNaN'),'*','Color',gre,'MarkerSize',10,'LineWidth',2);
% 
% hold off;
% ylim([0,1])
% ylabel('PLV','FontWeight','bold')
% xlabel('Frequency(Hz)','FontWeight','bold')


subplot(2,2,1)
axes('Units', 'Normalized', 'Position',[.33 .8 .12 .1])
box on
hold on
plot(mean(TTS_Pre_t_all,1, 'omitNaN'), mean(TTS_Pre_T_env_all,1, 'omitNaN'),'Color',blck, 'LineWidth',2);
plot(mean(TTS_Post_t_all,1, 'omitNaN'), mean(TTS_Post_T_env_all,1, 'omitNaN'),'Color',rd, 'LineWidth',2);
xlim([0.3,.35]);
ylim([-1.2,1.2]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

subplot(2,2,3)
axes('Units', 'Normalized', 'Position',[.33 .325 .12 .1])
box on
hold on
plot(mean(CA_Pre_t_all,1, 'omitNaN'), mean(CA_Pre_T_env_all,1, 'omitNaN'),'Color',blck, 'LineWidth',2);
plot(mean(CA_Post_t_all,1, 'omitNaN'), mean(CA_Post_T_env_all,1, 'omitNaN'),'Color',blu, 'LineWidth',2);
xlim([0.3,.35]);
ylim([-1.2,1.2]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

subplot(2,2,2)
axes('Units', 'Normalized', 'Position',[.77 .8 .12 .1])
box on
hold on
plot(mean(PTS_Pre_t_all,1, 'omitNaN'), mean(PTS_Pre_T_env_all,1, 'omitNaN'),'Color',blck, 'LineWidth',2);
plot(mean(PTS_Post_t_all,1, 'omitNaN'), mean(PTS_Post_T_env_all,1, 'omitNaN'),'Color',yel, 'LineWidth',2);
xlim([0.3,.35]);
ylim([-1.2,1.2]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

% subplot(2,2,4)
% axes('Units', 'Normalized', 'Position',[.77 .325 .12 .1])
% box on
% hold on
% plot(mean(GE_Pre_t_all,1, 'omitNaN'), mean(GE_Pre_T_env_all,1, 'omitNaN'),'Color',yel, 'LineWidth',2);
% plot(mean(GE_Post_t_all,1, 'omitNaN'), mean(GE_Post_T_env_all,1, 'omitNaN'),'Color',yel, 'LineWidth',2);
% xlim([0.3,.4]);
% ylim([-2,2]);
% yticks([-1,0,1])
% xlabel('Time(s)','FontWeight','bold');
% ylabel('Amplitude \muV','FontWeight','bold')
% hold off

subplot(2,2,4)
hold on; 
for i = 1:TTS_count
    plot([0,.5], [TTS_pre_ratio(i,1), TTS_post_ratio(i,1)], 'o-', 'Color', rd, 'linew',4, 'MarkerSize', 8)
end
for i = 1:CA_count
    plot([1.5,2], [CA_pre_ratio(i,1), CA_post_ratio(i,1)], 'o-', 'Color', blu, 'linew',4, 'MarkerSize', 8)
end
for i = 1:PTS_count
    plot([3,3.5], [PTS_pre_ratio(i,1), PTS_post_ratio(i,1)], 'o-', 'Color', yel, 'linew',4, 'MarkerSize', 8)
end
xlim([-.5,4])
xticks([.25, 1.75, 3.25])
xticklabels({'Syn','IHC','Complex'})
ylabel('\Sigma{Harmonics 1:5}')
set(gca, 'FontSize', 20, 'FontWeight', 'Bold')
title('Pre-Post', 'FontSize', 24, 'FontWeight', 'Bold');
grid on


%%
% 
cd 'Figures'
print -dpng -r600 Pre-Post-RAM
cd .. 