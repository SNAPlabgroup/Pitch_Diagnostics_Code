%% Script to plot all data collected

clear;

chins = {'Q412', 'Q424', 'Q426', 'Q430', 'Q431', 'Q427', 'Q428', 'Q421', 'Q425', 'Q443', 'Q422', 'Q440', 'Q441'};
group = {'TTS', 'TTS', 'TTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA'};

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
        [TTS_Pre_power_all(TTS_count, :),TTS_Pre_deltapows_all(TTS_count,:), ...
            TTS_Pre_elicitors_all(TTS_count,:),TTS_Pre_thresh_all(TTS_count,:)] = get_data_memr(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [TTS_Post_power_all(TTS_count, :),TTS_Post_deltapows_all(TTS_count,:), ...
            TTS_Post_elicitors_all(TTS_count,:),TTS_Post_thresh_all(TTS_count,:)] = get_data_memr(prefix, subj, 'TTS_2wksPost');
        
    elseif strcmp(group{k}, 'PTS')
        PTS_count = PTS_count+1; %add new chin to the list
        PTS_chin{PTS_count} = chins{k};
        
        % first get baseline data
        [PTS_Pre_power_all(PTS_count, :),PTS_Pre_deltapows_all(PTS_count,:), ...
            PTS_Pre_elicitors_all(PTS_count,:),PTS_Pre_thresh_all(PTS_count,:)] = get_data_memr(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [PTS_Post_power_all(PTS_count, :),PTS_Post_deltapows_all(PTS_count,:), ...
            PTS_Post_elicitors_all(PTS_count,:),PTS_Post_thresh_all(PTS_count,:)] = get_data_memr(prefix, subj, 'PTS_2wksPost');
        
        
    elseif strcmp(group{k}, 'CA')
        CA_count = CA_count+1; %add new chin to the list
        CA_chin{CA_count} = chins{k};
        
        % first get baseline data
        [CA_Pre_power_all(CA_count, :),CA_Pre_deltapows_all(CA_count,:), ...
            CA_Pre_elicitors_all(CA_count,:),CA_Pre_thresh_all(CA_count,:)] = get_data_memr(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [CA_Post_power_all(CA_count, :),CA_Post_deltapows_all(CA_count,:), ...
            CA_Post_elicitors_all(CA_count,:),CA_Post_thresh_all(CA_count,:)] = get_data_memr(prefix, subj, 'CA_2wksPost');
        
        
    elseif strcmp(group{k}, 'GE')
        GE_count = GE_count+1; %add new chin to the list
        GE_chin{GE_count} = chins{k};
        
        % first get baseline data
        [GE_Pre_power_all(GE_count, :),GE_Pre_deltapows_all(GE_count,:), ...
            GE_Pre_elicitors_all(GE_count,:),GE_Pre_thresh_all(GE_count,:)] = get_data_memr(prefix, subj, 'Baseline');
        
        % Then get 2wks Post data
        [GE_Post_power_all(GE_count, :),GE_Post_deltapows_all(GE_count,:), ...
            GE_Post_elicitors_all(GE_count,:),GE_Post_thresh_all(GE_count,:)] = get_data_memr(prefix, subj, 'GE_2wksPost');
        
        
    end
    
end



PTS_Pre_deltapows_all(2,:) = interp1(PTS_Pre_elicitors_all(2,:), PTS_Pre_deltapows_all(2,:), PTS_Pre_elicitors_all(1,:)); 
PTS_Pre_elicitors_all(2,:) = PTS_Pre_elicitors_all(1,:); 
%% Plot Data

% blck = [0.25, 0.25, 0.25];
% rd = [216, 27, 96]./255; %TTS
% blu = [30, 136, 229]./255; %CA
% yel = [255, 193, 7]./255; %PTS
% gre = [115, 177, 117]./255; %GE
% 
% i_blck = [0.75, 0.75, 0.75];
% i_rd = [216, 27, 96, 75]./255; %TTS
% i_blu = [30, 136, 229, 75]./255; %CA
% i_yel = [255, 193, 7, 75]./255; %PTS
% i_gre = [115, 177, 117, 57]./255; %GE

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

figure; %Spectral Domain
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])
% TTS
subplot(2,2,1)
hold on;
tts_count = sum(strcmp(group,'TTS')); 
% for i = 1:tts_count
%     plot(TTS_Pre_elicitors_all(i,:),TTS_Pre_deltapows_all(i,:),'Color',i_blck,'linewidth',3);
%     plot(TTS_Post_elicitors_all(i,:),TTS_Post_deltapows_all(i,:),'Color',i_rd,'linewidth',3);
% end
plot(mean(TTS_Pre_elicitors_all,1, 'omitNaN'),mean(TTS_Pre_deltapows_all,1, 'omitNaN'),'-o','Color',blck,'linewidth',4);
plot(mean(TTS_Post_elicitors_all,1, 'omitNaN'),mean(TTS_Post_deltapows_all,1, 'omitNaN'),'-o','Color',rd,'linewidth',4);

hold off;
ylim([0,5])
ylabel('\Delta Power')
xlabel('Elicitor Level (dB FPL)')
xlim([40,110])
set(gca, 'FontSize', 20)
title('Synaptopathy', 'FontSize', 24, 'Color', rd);
grid on; 

% Carbo
subplot(2,2,3)
hold on;
ca_count = sum(strcmp(group,'CA')); 
% for i = 1:ca_count
%     
%     plot(CA_Pre_elicitors_all(i,:),CA_Pre_deltapows_all(i,:),'Color',i_blck,'linewidth',3);
%     plot(CA_Post_elicitors_all(i,:),CA_Post_deltapows_all(i,:),'Color',i_blu,'linewidth',3);
%     
% end
plot(mean(CA_Pre_elicitors_all,1, 'omitNaN'),mean(CA_Pre_deltapows_all,1, 'omitNaN'),'-o','Color',blck,'linewidth',4);
plot(mean(CA_Post_elicitors_all,1, 'omitNaN'),mean(CA_Post_deltapows_all,1, 'omitNaN'),'-o','Color',blu,'linewidth',4);

hold off;
ylim([0,5])
xlim([40,110])
ylabel('\Delta Power')
xlabel('Elicitor Level (dB FPL)')
set(gca, 'FontSize', 20)
title('IHC Dysfunction', 'FontSize', 24, 'Color', blu);
grid on; 

% PTS
subplot(2,2,2)
hold on;
pts_count = sum(strcmp(group,'PTS')); 
% for i = 1:pts_count
%     
%     plot(PTS_Pre_elicitors_all(i,:),PTS_Pre_deltapows_all(i,:),'Color',i_blck,'linewidth',3);
%     plot(PTS_Post_elicitors_all(i,:),PTS_Post_deltapows_all(i,:),'Color',i_yel,'linewidth',3);
%     
% end
plot(mean(PTS_Pre_elicitors_all,1, 'omitNaN'),mean(PTS_Pre_deltapows_all,1, 'omitNaN'),'-o','Color',blck,'linewidth',4);

hold on;
plot(mean(PTS_Post_elicitors_all,1, 'omitNaN'),mean(PTS_Post_deltapows_all,1, 'omitNaN'),'-o','Color',yel,'linewidth',4);

hold off;
ylim([0,5])
ylabel('\Delta Power')
xlabel('Elicitor Level (dB FPL)')
xlim([40,110])
set(gca, 'FontSize', 20)
title('Complex SNHL', 'FontSize', 24, 'Color', yel);
% subplot(2,2,4)
% hold on;
% title('GE | RAM - 25% Duty Cycle','FontSize',14);
% plot(mean(GE_Pre_elicitors_all,1, 'omitNaN'),mean(GE_Pre_thresh_env_all,1, 'omitNaN'),'Color',blck,'linewidth',1.5);
% 
% hold on;
% plot(mean(GE_Post_elicitors_all,1, 'omitNaN'),mean(GE_Post_thresh_env_all,1, 'omitNaN'),'Color',gre,'linewidth',1.5);
% 
% hold off;
% ylim([0,5])
% ylabel('PLV','FontWeight','bold')
% xlabel('Frequency(Hz)','FontWeight','bold')


subplot(2,2,4)
hold on; 
for i = 1:tts_count
    plot([0,.5], [TTS_Pre_thresh_all(i), TTS_Post_thresh_all(i)], 'o-', 'Color', rd, 'linew',4)
end
for i = 1:ca_count
    plot([1.5,2], [CA_Pre_thresh_all(i), CA_Post_thresh_all(i)], 'o-', 'Color', blu, 'linew',4)
end
for i = 1:pts_count
    plot([3,3.5], [PTS_Pre_thresh_all(i), PTS_Post_thresh_all(i)], 'o-', 'Color', yel, 'linew',4)
end
xlim([-.5,4])
xticks([.25, 1.75, 3.25])
xticklabels({'Syn','IHC','Complex'})
ylabel('Threshold (dB FPL)')
set(gca, 'FontSize', 20)
title('Thresholds', 'FontSize', 24)
grid on;


%% Save Figures
cd 'Figures'
print -dpng -r600 Pre-Post-MEMR
cd .. 