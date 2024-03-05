%% Awake vs Sedated

clear;
% 
chins = {'Q453'}; %'Q431', 'Q424', 'Q426', 'Q428', 'Q421', 'Q441', 'Q422'}; %Q442 422 427 425 440 423

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
       
    end
    
    for awake = 0       %sedated
        
        s_count = s_count+1; %add new chin to the list
        s_chin{s_count} = chins{k};
        
        [s_Pre_oae(s_count, :),s_Pre_nf(s_count,:), ...
            s_Pre_f2(s_count,:),s_Pre_oaesum(s_count,:), ...
            s_Pre_centerfreq(s_count,:)] = get_data_DP(prefix, subj, 'Baseline', awake);

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
title('Swept DPOAEs','FontSize',16);

for i = 1:size(a_Pre_f2,1)
    plot(a_Pre_f2(i,:),a_Pre_oae(i,:), 'Color',pre_sed_ind,'linewidth',2);
    plot(a_Pre_f2(i,:),a_Pre_nf(i,:), '--', 'Color',pre_sed_ind,'linewidth',2);
    plot(a_Pre_centerfreq(i,:),a_Pre_oaesum(i,:),'o','Color',pre_sed_group,'MarkerSize',10);
end


for i = 1:size(s_Pre_f2,1)
    plot(s_Pre_f2(i,:),s_Pre_oae(i,:),'Color',post_sed_ind,'linewidth',2);
        plot(s_Pre_f2(i,:),s_Pre_nf(i,:),'--', 'Color',post_sed_ind,'linewidth',2);
    plot(s_Pre_centerfreq(i,:),s_Pre_oaesum(i,:),'x','Color',post_sed_ind,'MarkerSize',10,'LineWidth',2);
end


plot(mean(a_Pre_f2,1, 'omitNaN'),mean(a_Pre_oae,1, 'omitNaN'),'Color',pre_awake_group,'linewidth',5);
plot(mean(a_Pre_f2,1, 'omitNaN'),mean(a_Pre_nf,1, 'omitNaN'),'--', 'Color',pre_awake_group,'linewidth',1.5, 'HandleVisibility','off');
plot(mean(a_Pre_centerfreq,1, 'omitNaN'),mean(a_Pre_oaesum,1, 'omitNaN'),'*','Color',pre_awake_group,'MarkerSize',10,'LineWidth',2, 'HandleVisibility','off');

plot(mean(s_Pre_f2,1, 'omitNaN'),mean(s_Pre_oae,1, 'omitNaN'),'Color',post_awake_group,'linewidth',5);
plot(mean(s_Pre_f2,1, 'omitNaN'),mean(s_Pre_nf,1, 'omitNaN'),'--', 'Color',post_awake_group,'linewidth',1.5, 'HandleVisibility','off');
plot(mean(s_Pre_centerfreq,1, 'omitNaN'),mean(s_Pre_oaesum,1, 'omitNaN'),'*','Color',post_awake_group,'MarkerSize',10,'LineWidth',2, 'HandleVisibility','off');

hold off;
ylim([-40, 60])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
legend('Awake','Sedated')
set(gca, 'XScale', 'log', 'FontSize', 14)


%% Difference Plot

diff_oae = s_Pre_oae - a_Pre_oae; 
diff_sum = s_Pre_oaesum - a_Pre_oaesum; 

figure; %Spectral Domain
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])

hold on;
text(.6,-30, sprintf('n = %d', size(a_Pre_f2, 1) - sum(isnan(a_Pre_f2(:,1)))), 'FontSize', 12)
title('Difference w/ Sedation','FontSize',16);

for i = 1:size(a_Pre_f2,1)
    plot(a_Pre_f2(i,:),diff_oae(i,:), 'Color',pre_sed_ind,'linewidth',2);
    %plot(a_Pre_f2(i,:),a_Pre_nf(i,:), '--', 'Color',pre_sed_ind,'linewidth',2);
    plot(a_Pre_centerfreq(i,:),diff_sum(i,:),'o','Color',pre_sed_group,'MarkerSize',10);
end


plot(mean(a_Pre_f2,1, 'omitNaN'),mean(diff_oae,1, 'omitNaN'),'Color',pre_awake_group,'linewidth',5);
plot(mean(a_Pre_centerfreq,1, 'omitNaN'),mean(diff_sum,1, 'omitNaN'),'*','Color',pre_awake_group,'MarkerSize',10,'LineWidth',2, 'HandleVisibility','off');
plot(a_Pre_centerfreq(1,:), zeros(size(a_Pre_centerfreq(1,:))), '--', 'Color', pre_awake_group, 'linew', 3); 

hold off;
ylim([-30, 30])
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
legend('Awake','Sedated')
set(gca, 'XScale', 'log', 'FontSize', 14)


