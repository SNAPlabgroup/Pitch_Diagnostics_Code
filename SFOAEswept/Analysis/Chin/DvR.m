%% D v R

clear;

chins = {'Q412', 'Q423', 'Q424', 'Q426', 'Q430', 'Q431', 'Q427', 'Q428'};
group = {'TTS', 'TTS', 'TTS', 'TTS', 'CA', 'CA', 'PTS', 'PTS'};

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

%% Get Data

%% GET DATA
for k = 1:length(chins)
    
    subj = chins{k};
    
    if strcmp(group(k), 'TTS')
        TTS_count = TTS_count+1; %add new chin to the list
        TTS_chin{TTS_count} = chins{k};
        
        
        % first get baseline data
        [~, ~, ~,TTS_Pre_oaesum_DP(TTS_count,:), TTS_Pre_centerfreq_DP(TTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        [~, ~, ~,TTS_Pre_oaesum_SF(TTS_count,:), TTS_Pre_centerfreq_SF(TTS_count,:)] = get_data_SF(prefix, subj, 'Baseline', awake);
        
        % Then get 1day Post data
        [~, ~, ~,TTS_Post_1d_oaesum_DP(TTS_count,:), TTS_Post_1d_centerfreq_DP(TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_1dayPost', awake);
        [~, ~, ~,TTS_Post_1d_oaesum_SF(TTS_count,:), TTS_Post_1d_centerfreq_SF(TTS_count,:)] = get_data_SF(prefix, subj, 'TTS_1dayPost', awake);
        
        
        % Then get 2wks Post data
        [~, ~, ~,TTS_Post_oaesum_DP(TTS_count,:), TTS_Post_centerfreq_DP(TTS_count,:)] = get_data_DP(prefix, subj, 'TTS_2wksPost', awake);
        [~, ~, ~,TTS_Post_oaesum_SF(TTS_count,:), TTS_Post_centerfreq_SF(TTS_count,:)] = get_data_SF(prefix, subj, 'TTS_2wksPost', awake);
        
        
    elseif strcmp(group{k}, 'PTS')
        PTS_count = PTS_count+1; %add new chin to the list
        PTS_chin{PTS_count} = chins{k};
        
        % first get baseline data
        [~, ~, ~,PTS_Pre_oaesum_DP(PTS_count,:), PTS_Pre_centerfreq_DP(PTS_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        [~, ~, ~,PTS_Pre_oaesum_SF(PTS_count,:), PTS_Pre_centerfreq_SF(PTS_count,:)] = get_data_SF(prefix, subj, 'Baseline', awake);
        
        
        % Then get 2wks Post data
        [~, ~, ~,PTS_Post_oaesum_DP(PTS_count,:), PTS_Post_centerfreq_DP(PTS_count,:)] = get_data_DP(prefix, subj, 'PTS_2wksPost', awake);
        [~, ~, ~,PTS_Post_oaesum_SF(PTS_count,:), PTS_Post_centerfreq_SF(PTS_count,:)] = get_data_SF(prefix, subj, 'PTS_2wksPost', awake);
        
        
    elseif strcmp(group{k}, 'CA')
        CA_count = CA_count+1; %add new chin to the list
        CA_chin{CA_count} = chins{k};
        
        % first get baseline data
        [~, ~, ~,CA_Pre_oaesum_DP(CA_count,:), CA_Pre_centerfreq_DP(CA_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        [~, ~, ~,CA_Pre_oaesum_SF(CA_count,:), CA_Pre_centerfreq_SF(CA_count,:)] = get_data_SF(prefix, subj, 'Baseline', awake);
        
        
        % Then get 2wks Post data
        [~, ~, ~,CA_Post_oaesum_DP(CA_count,:), CA_Post_centerfreq_DP(CA_count,:)] = get_data_DP(prefix, subj, 'CA_2wksPost', awake);
        [~, ~, ~,CA_Post_oaesum_SF(CA_count,:), CA_Post_centerfreq_SF(CA_count,:)] = get_data_SF(prefix, subj, 'CA_2wksPost', awake);
        
        
        
    elseif strcmp(group{k}, 'GE')
        GE_count = GE_count+1; %add new chin to the list
        GE_chin{GE_count} = chins{k};
        
        % first get baseline data
        [~, ~, ~,GE_Pre_oaesum_DP(GE_count,:), GE_Pre_centerfreq_DP(GE_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        [~, ~, ~,GE_Pre_oaesum_SF(GE_count,:), GE_Pre_centerfreq_SF(GE_count,:)] = get_data_SF(prefix, subj, 'Baseline', awake);
        
        
        % Then get 2wks Post data
        [~, ~, ~,GE_Post_oaesum_DP(GE_count,:), GE_Post_centerfreq_DP(GE_count,:)] = get_data_DP(prefix, subj, 'GE_2wksPost', awake);
        [~, ~, ~,GE_Post_oaesum_SF(GE_count,:), GE_Post_centerfreq_SF(GE_count,:)] = get_data_SF(prefix, subj, 'GE_2wksPost', awake);
        
        
        
    end
    
end

%% Plot Data

blck = [0.25, 0.25, 0.25];
rd = [216, 27, 96]./255; %TTS
blu = [30, 136, 229]./255; %CA
yel = [255, 193, 7]./255; %PTS
gre = [115, 177, 117]./255; %GE


figure; %Spectral Domain
hold on; 
set(gcf, 'Units', 'inches', 'Position', [1, 1, 8, 7.5])
title('D vs R','FontSize',16);
xlim([-10,50])
ylim([-10, 50])
plot([-10, 50], [-10, 50], '--', 'color', [.5, .5 .5])
plot(TTS_Pre_oaesum_DP, TTS_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', rd, 'MarkerSize', 10)
plot(TTS_Post_1d_oaesum_DP, TTS_Post_1d_oaesum_SF, 'o', 'MarkerFaceColor', rd, 'MarkerEdgeColor', rd, 'MarkerSize', 10)
plot(TTS_Post_oaesum_DP, TTS_Post_oaesum_SF, '*', 'MarkerFaceColor', rd, 'MarkerEdgeColor', rd, 'MarkerSize', 10)

plot(PTS_Pre_oaesum_DP, PTS_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', yel, 'MarkerSize', 10)
plot(PTS_Post_oaesum_DP, PTS_Post_oaesum_SF, 'o', 'MarkerFaceColor', yel, 'MarkerEdgeColor', yel, 'MarkerSize', 10)

plot(CA_Pre_oaesum_DP, CA_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', blu, 'MarkerSize', 10)
plot(CA_Post_oaesum_DP, CA_Post_oaesum_SF, 'o', 'MarkerFaceColor', blu, 'MarkerEdgeColor', blu, 'MarkerSize', 10)

% plot(GE_Pre_oaesum_DP, GE_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', gre, 'MarkerSize', 10)
% plot(GE_Post_oaesum_DP, GE_Post_oaesum_SF, 'o', 'MarkerFaceColor', gre, 'MarkerEdgeColor', gre, 'MarkerSize', 10)

ylabel('SF Amplitude (dB EPL)','FontWeight','bold')
xlabel('DF Amplitude (dB EPL)','FontWeight','bold')
set(gca, 'FontSize', 14)


%% Separated
figure; %Spectral Domain
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])
% TTS
subplot(2,2,1)
hold on;
title('TTS','FontSize',16);
text(.71,45, sprintf('n = %d', size(TTS_Pre_oaesum_DP, 1) - sum(isnan(TTS_Pre_oaesum_DP(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(TTS_Post_oaesum_DP, 1) - sum(isnan(TTS_Post_oaesum_DP(:,1)))), 'Color', rd, 'FontSize', 12)
text(.71,35, sprintf('n = %d', size(TTS_Post_1d_oaesum_DP, 1) - sum(isnan(TTS_Post_1d_oaesum_DP(:,1)))), 'Color', rd, 'FontSize', 12)
hold on; 
xlim([-10,50])
ylim([-10, 50])
plot([-10, 50], [-10, 50], '--', 'color', [.5, .5 .5])
plot(TTS_Pre_oaesum_DP, TTS_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', rd, 'MarkerSize', 10)
plot(TTS_Post_1d_oaesum_DP, TTS_Post_1d_oaesum_SF, 'o', 'MarkerFaceColor', rd, 'MarkerEdgeColor', rd, 'MarkerSize', 10)
plot(TTS_Post_oaesum_DP, TTS_Post_oaesum_SF, '*', 'MarkerFaceColor', rd, 'MarkerEdgeColor', rd, 'MarkerSize', 10)
ylabel('SF Amplitude (dB EPL)','FontWeight','bold')
xlabel('DP Amplitude (dB EPL)','FontWeight','bold')
set(gca, 'FontSize', 14)

subplot(2,2,2)
hold on;
title('PTS','FontSize',16);
xlim([-10,50])
ylim([-10, 50])
plot([-10, 50], [-10, 50], '--', 'color', [.5, .5 .5])
text(.71,45, sprintf('n = %d', size(PTS_Pre_oaesum_DP, 1) - sum(isnan(PTS_Pre_oaesum_DP(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(PTS_Post_oaesum_DP, 1) - sum(isnan(PTS_Post_oaesum_DP(:,1)))), 'Color', yel, 'FontSize', 12)
plot(PTS_Pre_oaesum_DP, PTS_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', yel, 'MarkerSize', 10)
plot(PTS_Post_oaesum_DP, PTS_Post_oaesum_SF, 'o', 'MarkerFaceColor', yel, 'MarkerEdgeColor', yel, 'MarkerSize', 10)
ylabel('SF Amplitude (dB EPL)','FontWeight','bold')
xlabel('DP Amplitude (dB EPL)','FontWeight','bold')
set(gca, 'FontSize', 14)

subplot(2,2,3)
hold on;
title('CA','FontSize',16);
xlim([-10,50])
ylim([-10, 50])
plot([-10, 50], [-10, 50], '--', 'color', [.5, .5 .5])
text(.71,45, sprintf('n = %d', size(CA_Pre_oaesum_DP, 1) - sum(isnan(CA_Pre_oaesum_DP(:,1)))), 'FontSize', 12)
text(.71,40, sprintf('n = %d', size(CA_Post_oaesum_DP, 1) - sum(isnan(CA_Post_oaesum_DP(:,1)))), 'Color', blu, 'FontSize', 12)
plot(CA_Pre_oaesum_DP, CA_Pre_oaesum_SF, 'o', 'MarkerFaceColor', blck, 'MarkerEdgeColor', blu, 'MarkerSize', 10)
plot(CA_Post_oaesum_DP, CA_Post_oaesum_SF, 'o', 'MarkerFaceColor', blu, 'MarkerEdgeColor', blu, 'MarkerSize', 10)
ylabel('SF Amplitude (dB EPL)','FontWeight','bold')
xlabel('DP Amplitude (dB EPL)','FontWeight','bold')
set(gca, 'FontSize', 14)
