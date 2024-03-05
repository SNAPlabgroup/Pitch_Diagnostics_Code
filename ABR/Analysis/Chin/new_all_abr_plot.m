%% Plot ABR thresholds from new analysis
clear;
% chins = {'Q412', 'Q424', 'Q426', 'Q430', 'Q431', 'Q427', 'Q428', 'Q421', 'Q425', 'Q443', 'Q422', 'Q440', 'Q441'};
% group = {'TTS', 'TTS', 'TTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA', 'PTS', 'PTS', 'CA', 'CA'};

chin = []; % placeholder in case we want to only plot one chins data
thresholds = [1 2 3 4 5]; % placeholder while not loading real data

% Directories for individual user // currently for SH
datadir = '/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/ABR/Chin/';
cd([datadir, filesep, 'Baseline'])

if isempty(chin)
    all_chins = dir('Q*');
else
    all_chins = dir(chin);
end

% Get list of chins for each exposure
cd(datadir);
cd('CA_2wksPost')
ca_chins = dir('Q*');
ca_chins = {ca_chins.name};
cd(datadir);
cd('PTS_2wksPost')
pts_chins = dir('Q*');
pts_chins = {pts_chins.name};
cd(datadir);
cd('TTS_2wksPost')
tts_chins = dir('Q*');
tts_chins = {tts_chins.name};
cd(datadir);
cd('GE_2wksPost')
ge_chins = dir('Q*');
ge_chins = {ge_chins.name};

cd(datadir)

% Initialize results/data matrices
freq = [.5, 1, 2, 4, 8]; % kHz
baseline = zeros(numel(all_chins),numel(freq));
post = zeros(numel(all_chins),numel(freq));
exp = [];

for k = 1:numel(all_chins)
    chin = all_chins(k).name;
    
    % Get Baseline data
    cd([all_chins(k).folder,filesep, chin, filesep, 'Processed'])
    % load(XXX);
    baseline(k,:) = 10*rand(1,5); %thresholds;
    
    cd(datadir);
    emptyFlag = 0;
    if sum(strcmp(chin, tts_chins)>0)
        cd(fullfile('TTS_2wksPost', chin, 'Processed'))
        exp{k,1} = 'TTS';
    elseif sum(strcmp(chin, pts_chins)>0)
        cd(fullfile('PTS_2wksPost', chin, 'Processed'))
        exp{k,1} = 'PTS';
    elseif sum(strcmp(chin, ca_chins)>0)
        cd(fullfile('CA_2wksPost', chin, 'Processed'))
        exp{k,1} = 'CA';
    elseif sum(strcmp(chin, ge_chins)>0)
        cd(fullfile('GE_2wksPost', chin, 'Processed'))
        exp{k,1} = 'GE';
    else
        exp{k,1} = 'NA';
        emptyFlag = 1;
    end
    
    if ~emptyFlag
        %load(XXX)
        post(k,:) = 20*rand(1,5); %thresholds;
    end
    
end

%% Plot Data

blck = [0.25, 0.25, 0.25];
rd = [194 106 119]./255; %TTS
blu = [148 203 236]./255; %CA
yel = [220 205 125]./255; %PTS
gre = [93 168 153]./255; %GE

i_blck = [0.25, 0.25, .25, 75];
i_rd = [194 106 119 75]./255; %TTS
i_blu = [148 203 236 75]./255; %CA
i_yel = [220 205 125 75]./255; %PTS
i_gre = [93 168 153 57]./255; %GE

i_cols = [i_blck; i_rd; i_blu; i_yel; i_gre]; 
cols = [blck; rd; blu; yel; gre]; 
groups = {'NH', 'TTS', 'CA', 'PTS', 'GE'}; 
subp = [0 1 3 2 4]'; 

figure;
hold on; 
set(gcf, 'Units', 'inches', 'Position', [1, 1, 16, 12])

if size(all_chins,1) == 1
    plot(freq,baseline, '-o','Color',blck,'linewidth',4, 'MarkerSize', 8)
    plot(freq,baseline, '-o','Color',rd,'linewidth',4, 'MarkerSize', 8)
else
    for j = 1:numel(all_chins)
        grp = strcmp(exp{j}, groups);
        if sum(grp) > 0
            subplot(2,2,grp * subp)
            hold on;
            plot(freq,baseline(j,:), '-o','Color',i_blck(1,1:3),'linewidth',4, 'MarkerSize', 8)
            plot(freq,post(j,:), '-o','Color',grp * i_cols,'linewidth',4, 'MarkerSize', 8)
           text(9,baseline(j,5), all_chins(j).name, 'Units', 'Data', 'Color', i_blck(1,1:3))
           text(9,post(j,5), all_chins(j).name, 'Units', 'Data', 'Color',grp * i_cols)
        end
    end
end
    % TTS
    subplot(2,2,1)
    hold on;
    title('Synaptopathy');
    
    % for i = 1:tts_count
    %     plot(freq,TTS_pre(i,:),'Color',i_blck,'linewidth',3);
    %     plot(freq,TTS_post(i,:),'Color',i_rd,'linewidth',3);
    % end
%     plot(freq,mean(TTS_pre,1, 'omitNaN'),'-o','Color',blck,'linewidth',4, 'MarkerSize', 8);
%     plot(freq,mean(TTS_post,1, 'omitNaN'),'-o','Color',rd,'linewidth',4, 'MarkerSize', 8);
%     
    hold off;
    ylim([0,60])
    ylabel('Threshold (dB SPL)')
    xlabel('Frequency (kHz)')
    xlim([0.5, 8])
    xticks(freq)
    set(gca, 'FontSize', 20, 'XScale', 'log')
    title('Synaptopathy', 'FontSize', 24, 'Color', rd);
    grid on
    
    % Carbo
    subplot(2,2,3)
    hold on;
    title('IHC Dysfunction');
    % for i = 1:ca_count
    %     plot(freq,CA_pre(i,:),'Color',i_blck,'linewidth',3);
    %     plot(freq,CA_post(i,:),'Color',i_blu,'linewidth',3);
    % end
%     plot(freq,mean(CA_pre,1, 'omitNaN'),'-o','Color',blck,'linewidth',4, 'MarkerSize', 8);
%     plot(freq,mean(CA_post,1, 'omitNaN'),'-o','Color',blu,'linewidth',4, 'MarkerSize', 8);
%     
    hold off;
    ylim([0,60])
    ylabel('Threshold (dB SPL)')
    xlabel('Frequency (kHz)')
    xlim([0.5, 8])
    xticks(freq)
    set(gca, 'FontSize', 20, 'XScale', 'log')
    title('IHC Dysfunction', 'FontSize', 24, 'Color', blu);
    grid on
    
    % PTS
    subplot(2,2,2)
    hold on;
    % for i = 1:pts_count
    %     plot(freq,PTS_pre(i,:),'Color',i_blck,'linewidth',3);
    %     plot(freq,PTS_post(i,:),'Color',i_yel,'linewidth',3);
    % end
%     plot(freq,mean(PTS_pre,1, 'omitNaN'),'-o','Color',blck,'linewidth',4, 'MarkerSize', 8);
%     plot(freq,mean(PTS_post,1, 'omitNaN'),'-o','Color',yel,'linewidth',4, 'MarkerSize', 8);
%     
    hold off;
    ylim([0,60])
    ylabel('Threshold (dB SPL)')
    xlabel('Frequency (kHz)')
    xlim([0.5, 8])
    xticks(freq)
    set(gca, 'FontSize', 20, 'XScale', 'log')
    title('Complex SNHL', 'FontSize', 24, 'Color', yel);
    grid on;
    
    subplot(2,2,4)
    hold on;
%     for i = 1:tts_count
%         plot([0,.5], [mean(TTS_pre(i,:)), mean(TTS_post(i,:))], 'o-', 'Color', rd, 'linew',4, 'MarkerSize', 8)
%     end
%     for i = 1:ca_count
%         plot([1.5,2], [mean(CA_pre(i,:)), mean(CA_post(i,:))], 'o-', 'Color', blu, 'linew',4, 'MarkerSize', 8)
%     end
%     for i = 1:pts_count
%         plot([3,3.5], [mean(PTS_pre(i,:)), mean(PTS_post(i,:))], 'o-', 'Color', yel, 'linew',4, 'MarkerSize', 8)
%     end
    ylim([0,60])
    ylabel('Threshold (dB SPL)')
    xlabel('Frequency (kHz)')
    xlim([0.5, 8])
    xticks(freq)
    set(gca, 'FontSize', 20, 'XScale', 'log')
    title('Gentamicin', 'FontSize', 24, 'Color', gre);
    grid on;
    
    

% %%
% cd 'Figures'
% print -dpng -r600 Pre-Post-ABR
% cd ..