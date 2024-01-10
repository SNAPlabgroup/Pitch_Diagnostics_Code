%% Awake vs Sedated

clear;
% 
chins = {'Q431'};
group = {'CA'};
% 
chins = {'Q424', 'Q426', 'Q412', 'Q423'};
group = {'TTS'};

% chins = {'Q427', 'Q428'};
% group = {'PTS'};


location = 0; % 0 == mac, 1 == Desktop, 2 == SNAPlab

uname = 'samhauser';

if location == 1 % School
    prefix = 'F:\';
elseif location == 0 % Mac
    prefix = ['/Volumes/SNH/'];
end

a_count = 0;

s_count = 0;

codedir = pwd;

%% GET DATA
for k = 1:length(chins)
    
    subj = chins{k};
    
    for awake = 1       % awake
        
        a_count = a_count+1; %add new chin to the list
        a_chin{a_count} = chins{k};
        
        % first get baseline data
        [a_Pre_oae(a_count, :),a_Pre_nf(a_count,:), ...
            a_Pre_f2(a_count,:),a_Pre_oaesum(a_count,:), ...
            a_Pre_centerfreq(a_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 2wks Post data
        [a_Post_oae(a_count, :),a_Post_nf(a_count,:), ...
            a_Post_f2(a_count,:),a_Post_oaesum(a_count,:), ...
            a_Post_centerfreq(a_count,:)] = get_data_DP(prefix, subj, [group{1},'_2wksPost'], awake);
       
    end
    
    for awake = 0       %sedated
        
        s_count = s_count+1; %add new chin to the list
        s_chin{s_count} = chins{k};
        
        % first get baseline data
        [s_Pre_oae(s_count, :),s_Pre_nf(s_count,:), ...
            s_Pre_f2(s_count,:),s_Pre_oaesum(s_count,:), ...
            s_Pre_centerfreq(s_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);
        
        % Then get 2wks Post data
        [s_Post_oae(s_count, :),s_Post_nf(s_count,:), ...
            s_Post_f2(s_count,:),s_Post_oaesum(s_count,:), ...
            s_Post_centerfreq(s_count,:)] = get_data_DP(prefix, subj, [group{1},'_2wksPost'], awake);
        
    end
end


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

hold on;
text(.6,-30, sprintf('n = %d', size(a_Pre_f2, 1) - sum(isnan(a_Pre_f2(:,1)))), 'FontSize', 12)
text(.6,-32.5, sprintf('n = %d', size(a_Post_f2, 1) - sum(isnan(a_Post_f2(:,1)))), 'Color', post_awake_group, 'FontSize', 12)
title(sprintf('%s | Swept DPOAEs', group{1}),'FontSize',16);
% 
% for i = 1:size(a_Pre_f2,1)
%     plot(a_Pre_f2(i,:),a_Pre_oae(i,:), 'Color',pre_awake_ind,'linewidth',2.5);
%     plot(a_Pre_centerfreq(i,:),a_Pre_oaesum(i,:),'o','Color',pre_awake_group,'MarkerSize',10);
% end
% 
% for i = 1:size(a_Post_f2,1)
%     plot(a_Post_f2(i,:),a_Pre_oae(i,:), 'Color',post_awake_ind,'linewidth',2.5);
%     plot(a_Post_centerfreq(i,:),a_Post_oaesum(i,:),'o','Color',post_awake_group,'MarkerSize',10,'LineWidth',2);
% end
% 
% for i = 1:size(s_Pre_f2,1)
%     plot(s_Pre_f2(i,:),s_Pre_oae(i,:),'Color',pre_sed_ind,'linewidth',2.5);
%     plot(s_Pre_centerfreq(i,:),s_Pre_oaesum(i,:),'*','Color',pre_sed_ind,'MarkerSize',10,'LineWidth',2);
% end
% 
% for i = 1:size(s_Pre_f2,1)
%     plot(s_Post_f2(i,:),s_Post_oae(i,:),'Color',post_sed_ind,'linewidth',2.5);
%     plot(s_Post_centerfreq(i,:),s_Post_oaesum(i,:),'o','Color',post_sed_group,'MarkerSize',10,'LineWidth',2);
% end

plot(mean(a_Pre_f2,1, 'omitNaN'),mean(a_Pre_oae,1, 'omitNaN'),'Color',pre_awake_group,'linewidth',5);
plot(mean(a_Pre_f2,1, 'omitNaN'),mean(a_Pre_nf,1, 'omitNaN'),'--', 'Color',pre_awake_group,'linewidth',1.5);
plot(mean(a_Pre_centerfreq,1, 'omitNaN'),mean(a_Pre_oaesum,1, 'omitNaN'),'*','Color',pre_awake_group,'MarkerSize',10,'LineWidth',2);
plot(mean(a_Post_f2,1, 'omitNaN'),mean(a_Post_oae,1, 'omitNaN'),'Color',post_awake_group,'linewidth',5);
plot(mean(a_Post_f2,1, 'omitNaN'),mean(a_Post_nf,1, 'omitNaN'),'--', 'Color',post_awake_group,'linewidth',1.5);
plot(mean(a_Post_centerfreq,1, 'omitNaN'),mean(a_Post_oaesum,1, 'omitNaN'),'*','Color',post_awake_group,'MarkerSize',10,'LineWidth',2);

plot(mean(s_Pre_f2,1, 'omitNaN'),mean(s_Pre_oae,1, 'omitNaN'),'Color',pre_sed_group,'linewidth',5);
plot(mean(s_Pre_f2,1, 'omitNaN'),mean(s_Pre_nf,1, 'omitNaN'),'--', 'Color',pre_sed_group,'linewidth',1.5);
plot(mean(s_Pre_centerfreq,1, 'omitNaN'),mean(s_Pre_oaesum,1, 'omitNaN'),'*','Color',pre_sed_group,'MarkerSize',10,'LineWidth',2);
plot(mean(s_Post_f2,1, 'omitNaN'),mean(s_Post_oae,1, 'omitNaN'),'Color',post_sed_group,'linewidth',5);
plot(mean(s_Post_f2,1, 'omitNaN'),mean(s_Post_nf,1, 'omitNaN'),'--', 'Color',post_sed_group,'linewidth',1.5);
plot(mean(s_Post_centerfreq,1, 'omitNaN'),mean(s_Post_oaesum,1, 'omitNaN'),'*','Color',post_sed_group,'MarkerSize',10,'LineWidth',2);

hold off;
ylim([-40, 50])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
set(gca, 'XScale', 'log', 'FontSize', 14)



