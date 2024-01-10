%% Awake vs Sedated

clear;

chins = {'Q424', 'Q426', 'Q412'};
group = {'TTS'}; 

location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

a_TTS_count = 0;
a_PTS_count = 0;
a_CA_count = 0;
a_GE_count = 0;

s_TTS_count = 0;
s_PTS_count = 0;
s_CA_count = 0;
s_GE_count = 0;

codedir = pwd;


%% GET DATA
for k = 1:length(chins)
    
    subj = chins{k};
    for awake = 1       % awake
        
        if strcmp(group(k), 'TTS')
            a_TTS_count = a_TTS_count+1; %add new chin to the list
            a_TTS_chin{a_TTS_count} = chins{k};
            
            
            % first get baseline data
            [a_TTS_Pre_oae(a_TTS_count, :),a_TTS_Pre_nf(a_TTS_count,:), ...
                a_TTS_Pre_f2(a_TTS_count,:),a_TTS_Pre_oaesum(a_TTS_count,:), ...
                a_TTS_Pre_centerfreq(a_TTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 1day Post data
            [a_TTS_Post_1d_oae(a_TTS_count, :),a_TTS_Post_1d_nf(a_TTS_count,:), ...
                a_TTS_Post_1d_f2(a_TTS_count,:),a_TTS_Post_1d_oaesum(a_TTS_count,:), ...
                a_TTS_Post_1d_centerfreq(a_TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_1dayPost', awake);
            
            % Then get 2wks Post data
            [a_TTS_Post_oae(a_TTS_count, :),a_TTS_Post_nf(a_TTS_count,:), ...
                a_TTS_Post_f2(a_TTS_count,:),a_TTS_Post_oaesum(a_TTS_count,:), ...
                a_TTS_Post_centerfreq(a_TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_2wksPost', awake);
            
        elseif strcmp(group{k}, 'PTS')
            a_PTS_count = a_PTS_count+1; %add new chin to the list
            a_PTS_chin{a_PTS_count} = chins{k};
            
            % first get baseline data
            [a_PTS_Pre_oae(a_PTS_count, :),a_PTS_Pre_nf(a_PTS_count,:), ...
                a_PTS_Pre_f2(a_PTS_count,:),a_PTS_Pre_oaesum(a_PTS_count,:), ...
                a_PTS_Pre_centerfreq(a_PTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [a_PTS_Post_oae(a_PTS_count, :),a_PTS_Post_nf(a_PTS_count,:), ...
                a_PTS_Post_f2(a_PTS_count,:),a_PTS_Post_oaesum(a_PTS_count,:), ...
                a_PTS_Post_centerfreq(a_PTS_count,:)] = get_data_DP(prefix, subj, 'PTS_2wksPost', awake);
            
            
        elseif strcmp(group{k}, 'CA')
            a_CA_count = a_CA_count+1; %add new chin to the list
            a_CA_chin{a_CA_count} = chins{k};
            
            % first get baseline data
            [a_CA_Pre_oae(a_CA_count, :),a_CA_Pre_nf(a_CA_count,:), ...
                a_CA_Pre_f2(a_CA_count,:),a_CA_Pre_oaesum(a_CA_count,:), ...
                a_CA_Pre_centerfreq(a_CA_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [a_CA_Post_oae(a_CA_count, :),a_CA_Post_nf(a_CA_count,:), ...
                a_CA_Post_f2(a_CA_count,:),a_CA_Post_oaesum(a_CA_count,:), ...
                a_CA_Post_centerfreq(a_CA_count,:)] = get_data_DP(prefix, subj, 'CA_2wksPost', awake);
            
            
        elseif strcmp(group{k}, 'GE')
            a_GE_count = a_GE_count+1; %add new chin to the list
            a_GE_chin{a_GE_count} = chins{k};
            
            % first get baseline data
            [a_GE_Pre_oae(a_GE_count, :),a_GE_Pre_nf(a_GE_count,:), ...
                a_GE_Pre_f2(a_GE_count,:),a_GE_Pre_oaesum(a_GE_count,:), ...
                a_GE_Pre_centerfreq(a_GE_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [a_GE_Post_oae(a_GE_count, :),a_GE_Post_nf(a_GE_count,:), ...
                a_GE_Post_f2(a_GE_count,:),a_GE_Post_oaesum(a_GE_count,:), ...
                a_GE_Post_centerfreq(a_GE_count,:)] = get_data_DP(prefix, subj, 'GE_2wksPost', awake);
            
            
        end
    end
    
    for awake = 0       %sedated
        
        if strcmp(group(k), 'TTS')
            s_TTS_count = s_TTS_count+1; %add new chin to the list
            s_TTS_chin{s_TTS_count} = chins{k};
            
            
            % first get baseline data
            [s_TTS_Pre_oae(s_TTS_count, :),s_TTS_Pre_nf(s_TTS_count,:), ...
                s_TTS_Pre_f2(s_TTS_count,:),s_TTS_Pre_oaesum(s_TTS_count,:), ...
                s_TTS_Pre_centerfreq(s_TTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 1day Post data
            [s_TTS_Post_1d_oae(s_TTS_count, :),s_TTS_Post_1d_nf(s_TTS_count,:), ...
                s_TTS_Post_1d_f2(s_TTS_count,:),s_TTS_Post_1d_oaesum(s_TTS_count,:), ...
                s_TTS_Post_1d_centerfreq(s_TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_1dayPost', awake);
            
            % Then get 2wks Post data
            [s_TTS_Post_oae(s_TTS_count, :),s_TTS_Post_nf(s_TTS_count,:), ...
                s_TTS_Post_f2(s_TTS_count,:),s_TTS_Post_oaesum(s_TTS_count,:), ...
                s_TTS_Post_centerfreq(s_TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_2wksPost', awake);
            
        elseif strcmp(group{k}, 'PTS')
            s_PTS_count = s_PTS_count+1; %add new chin to the list
            s_PTS_chin{s_PTS_count} = chins{k};
            
            % first get baseline data
            [s_PTS_Pre_oae(s_PTS_count, :),s_PTS_Pre_nf(s_PTS_count,:), ...
                s_PTS_Pre_f2(s_PTS_count,:),s_PTS_Pre_oaesum(s_PTS_count,:), ...
                s_PTS_Pre_centerfreq(s_PTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [s_PTS_Post_oae(s_PTS_count, :),s_PTS_Post_nf(s_PTS_count,:), ...
                s_PTS_Post_f2(s_PTS_count,:),s_PTS_Post_oaesum(s_PTS_count,:), ...
                s_PTS_Post_centerfreq(s_PTS_count,:)] = get_data_DP(prefix, subj, 'PTS_2wksPost', awake);
            
            
        elseif strcmp(group{k}, 'CA')
            a_CA_count = s_CA_count+1; %add new chin to the list
            s_CA_chin{s_CA_count} = chins{k};
            
            % first get baseline data
            [s_CA_Pre_oae(s_CA_count, :),s_CA_Pre_nf(s_CA_count,:), ...
                s_CA_Pre_f2(s_CA_count,:),s_CA_Pre_oaesum(s_CA_count,:), ...
                s_CA_Pre_centerfreq(s_CA_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [s_CA_Post_oae(s_CA_count, :),s_CA_Post_nf(s_CA_count,:), ...
                s_CA_Post_f2(s_CA_count,:),s_CA_Post_oaesum(s_CA_count,:), ...
                s_CA_Post_centerfreq(s_CA_count,:)] = get_data_DP(prefix, subj, 'CA_2wksPost', awake);
            
            
        elseif strcmp(group{k}, 'GE')
            s_GE_count = s_GE_count+1; %add new chin to the list
            s_GE_chin{s_GE_count} = chins{k};
            
            % first get baseline data
            [s_GE_Pre_oae(s_GE_count, :),s_GE_Pre_nf(s_GE_count,:), ...
                s_GE_Pre_f2(s_GE_count,:),s_GE_Pre_oaesum(s_GE_count,:), ...
                s_GE_Pre_centerfreq(s_GE_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
            
            % Then get 2wks Post data
            [s_GE_Post_oae(s_GE_count, :),s_GE_Post_nf(s_GE_count,:), ...
                s_GE_Post_f2(s_GE_count,:),s_GE_Post_oaesum(s_GE_count,:), ...
                s_GE_Post_centerfreq(s_GE_count,:)] = get_data_DP(prefix, subj, 'GE_2wksPost', awake);
            
            
        end
    end
    
    
    
end

%% Get averages

%% Plot Data

pre_awake_group= [77, 77, 77]./255; %pre awake
post_awake_group = [178, 24, 43]./255; %post awake
pre_sed_group = [186, 186, 186]./255; %pre sedated
post_sed_group = [244, 165, 130]./255; %post sedated

pre_awake_ind= [135, 135, 135]./255; %pre awake
post_awake_ind = [214, 96, 77]./255; %post awake
pre_sed_ind = [224, 224, 224]./255; %pre sedated
post_sed_ind = [253, 219, 199]./255; %post sedated

figure; %Spectral Domain
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])


% PTS
hold on;
text(.71,45, sprintf('n = %d', size(a_PTS_Pre_f2, 1) - sum(isnan(a_PTS_Pre_f2(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(a_PTS_Post_f2, 1) - sum(isnan(a_PTS_Post_f2(:,1)))), 'Color', post_awake_group, 'FontSize', 12)
title('PTS | Swept DPOAEs','FontSize',16);

plot(mean(a_PTS_Pre_f2,1, 'omitNaN'),mean(a_PTS_Pre_oae,1, 'omitNaN'),'Color',pre_awake_group,'linewidth',1.5);
plot(mean(a_PTS_Pre_f2,1, 'omitNaN'),mean(a_PTS_Pre_nf,1, 'omitNaN'),'--', 'Color',pre_awake_group,'linewidth',1.5);
plot(mean(a_PTS_Pre_centerfreq,1, 'omitNaN'),mean(a_PTS_Pre_oaesum,1, 'omitNaN'),'*','Color',pre_awake_group,'MarkerSize',10,'LineWidth',2);

for i = 1:size(a_PTS_Pre_f2,1)
    plot(a_PTS_Pre_f2(i,:),a_PTS_Pre_oae(i,:),':', 'Color',pre_awake_ind,'linewidth',1.5);
    plot(a_PTS_Pre_centerfreq(i,:),a_PTS_Pre_oaesum(i,:),'*','Color',pre_awake_group,'MarkerSize',10,'LineWidth',2);
end

hold on;
plot(mean(a_PTS_Post_f2,1, 'omitNaN'),mean(a_PTS_Post_oae,1, 'omitNaN'),'Color',post_awake_group,'linewidth',1.5);
plot(mean(a_PTS_Post_f2,1, 'omitNaN'),mean(a_PTS_Post_nf,1, 'omitNaN'),'--', 'Color',post_awake_group,'linewidth',1.5);
plot(mean(a_PTS_Post_centerfreq,1, 'omitNaN'),mean(a_PTS_Post_oaesum,1, 'omitNaN'),'*','Color',post_awake_group,'MarkerSize',10,'LineWidth',2);

for i = 1:size(a_PTS_Post_f2,1)
    plot(a_PTS_Post_f2(i,:),a_PTS_Pre_oae(i,:),':', 'Color',post_awake_ind,'linewidth',1.5);
    plot(a_PTS_Post_centerfreq(i,:),a_PTS_Post_oaesum(i,:),'*','Color',post_awake_group,'MarkerSize',10,'LineWidth',2);
end

hold on; 
plot(mean(s_PTS_Pre_f2,1, 'omitNaN'),mean(s_PTS_Pre_oae,1, 'omitNaN'),'Color',pre_sed_group,'linewidth',1.5);
plot(mean(s_PTS_Pre_f2,1, 'omitNaN'),mean(s_PTS_Pre_nf,1, 'omitNaN'),'--', 'Color',pre_sed_group,'linewidth',1.5);
plot(mean(s_PTS_Pre_centerfreq,1, 'omitNaN'),mean(s_PTS_Pre_oaesum,1, 'omitNaN'),'*','Color',pre_sed_group,'MarkerSize',10,'LineWidth',2);

for i = 1:size(s_PTS_Pre_f2,1)
    plot(s_PTS_Pre_f2(i,:),s_PTS_Pre_oae(i,:),':', 'Color',pre_sed_ind,'linewidth',1.5);
    plot(s_PTS_Pre_centerfreq(i,:),s_PTS_Pre_oaesum(i,:),'*','Color',pre_sed_ind,'MarkerSize',10,'LineWidth',2);
end

hold on;
plot(mean(s_PTS_Post_f2,1, 'omitNaN'),mean(s_PTS_Post_oae,1, 'omitNaN'),'Color',post_sed_group,'linewidth',1.5);
plot(mean(s_PTS_Post_f2,1, 'omitNaN'),mean(s_PTS_Post_nf,1, 'omitNaN'),'--', 'Color',post_sed_group,'linewidth',1.5);
plot(mean(s_PTS_Post_centerfreq,1, 'omitNaN'),mean(s_PTS_Post_oaesum,1, 'omitNaN'),'*','Color',post_sed_group,'MarkerSize',10,'LineWidth',2);

for i = 1:size(s_PTS_Pre_f2,1)
    plot(s_PTS_Post_f2(i,:),s_PTS_Post_oae(i,:),':', 'Color',post_sed_ind,'linewidth',1.5);
    plot(s_PTS_Post_centerfreq(i,:),s_PTS_Post_oaesum(i,:),'*','Color',post_sed_group,'MarkerSize',10,'LineWidth',2);
end

hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 14)




%% %%

% TTS
subplot(2,2,1)
hold on;
text(.71,45, sprintf('n = %d', size(TTS_Pre_f2, 1) - sum(isnan(TTS_Pre_f2(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(TTS_Post_f2, 1) - sum(isnan(TTS_Post_f2(:,1)))), 'Color', rd, 'FontSize', 12)
text(.71,35, sprintf('n = %d', size(TTS_Post_1d_f2, 1) - sum(isnan(TTS_Post_1d_f2(:,1)))), 'Color', rd, 'FontSize', 12)
title('TTS | Swept DPOAEs','FontSize',16);
plot(mean(TTS_Pre_f2,1, 'omitNaN'),mean(s_TTS_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',1.5);
plot(mean(TTS_Pre_f2,1, 'omitNaN'),mean(TTS_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',1.5);
plot(mean(TTS_Pre_centerfreq,1, 'omitNaN'),mean(TTS_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',2);

hold on;
plot(mean(TTS_Post_1d_f2,1, 'omitNaN'),mean(s_TTS_Post_1d_oae,1, 'omitNaN'),':', 'Color',rd,'linewidth',1.5);
plot(mean(TTS_Post_1d_f2,1, 'omitNaN'),mean(TTS_Post_1d_nf,1, 'omitNaN'),':', 'Color',rd,'linewidth',1.5);
plot(mean(TTS_Post_1d_centerfreq,1, 'omitNaN'),mean(TTS_Post_1d_oaesum,1, 'omitNaN'),'*','Color',rd,'MarkerSize',10,'LineWidth',2);

hold on;
plot(mean(TTS_Post_f2,1, 'omitNaN'),mean(TTS_Post_oae,1, 'omitNaN'),'Color',rd,'linewidth',1.5);
plot(mean(TTS_Post_f2,1, 'omitNaN'),mean(TTS_Post_nf,1, 'omitNaN'),'--', 'Color',rd,'linewidth',1.5);
plot(mean(TTS_Post_centerfreq,1, 'omitNaN'),mean(TTS_Post_oaesum,1, 'omitNaN'),'*','Color',rd,'MarkerSize',10,'LineWidth',2);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 14)

% CA
subplot(2,2,3)
hold on;
text(.71,45, sprintf('n = %d', size(CA_Pre_f2, 1) - sum(isnan(CA_Pre_f2(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(CA_Post_f2, 1) - sum(isnan(CA_Post_f2(:,1)))), 'Color', blu, 'FontSize', 12)
title('CA | Swept DPOAEs','FontSize',16);
plot(mean(CA_Pre_f2,1, 'omitNaN'),mean(CA_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',1.5);
plot(mean(CA_Pre_f2,1, 'omitNaN'),mean(CA_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',1.5);
plot(mean(CA_Pre_centerfreq,1, 'omitNaN'),mean(CA_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',2);

hold on;
plot(mean(CA_Post_f2,1, 'omitNaN'),mean(CA_Post_oae,1, 'omitNaN'),'Color',blu,'linewidth',1.5);
plot(mean(CA_Post_f2,1, 'omitNaN'),mean(CA_Post_nf,1, 'omitNaN'),'--', 'Color',blu,'linewidth',1.5);
plot(mean(CA_Post_centerfreq,1, 'omitNaN'),mean(CA_Post_oaesum,1, 'omitNaN'),'*','Color',blu,'MarkerSize',10,'LineWidth',2);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 14)

subplot(2,2,4)
hold on;
title('GE | Swept DPOAEs','FontSize',16);
text(.71,45, sprintf('n = %d', size(GE_Pre_f2, 1) - sum(isnan(GE_Pre_f2(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(GE_Post_f2, 1) - sum(isnan(GE_Post_f2(:,1)))), 'Color', gre, 'FontSize', 12)
plot(mean(GE_Pre_f2,1, 'omitNaN'),mean(a_GE_Pre_oae,1, 'omitNaN'),'Color',blck,'linewidth',1.5);
plot(mean(GE_Pre_f2,1, 'omitNaN'),mean(GE_Pre_nf,1, 'omitNaN'),'--', 'Color',blck,'linewidth',1.5);
plot(mean(GE_Pre_centerfreq,1, 'omitNaN'),mean(GE_Pre_oaesum,1, 'omitNaN'),'*','Color',blck,'MarkerSize',10,'LineWidth',2);

hold on;
plot(mean(GE_Post_f2,1, 'omitNaN'),mean(GE_Post_oae,1, 'omitNaN'),'Color',gre,'linewidth',1.5);
plot(mean(GE_Post_f2,1, 'omitNaN'),mean(GE_Post_nf,1, 'omitNaN'),'--', 'Color',gre,'linewidth',1.5);
plot(mean(GE_Post_centerfreq,1, 'omitNaN'),mean(GE_Post_oaesum,1, 'omitNaN'),'*','Color',gre,'MarkerSize',10,'LineWidth',2);
hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 14)


