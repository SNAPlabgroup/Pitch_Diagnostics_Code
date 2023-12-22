%% WBMEMR_Analysis

% DPOAE swept Analysis
% Author: Samantha Hauser modified code from Hari Bharadwaj
% Created: August 2023
% Last Updated: August 23, 2023
% Purpose:
% Helpful info:

%% Import data
cwd = pwd;
cd(datapath)

datafile = dir(fullfile(cd,['wbmemr*.mat']));
if length(datafile) < 1
    fprintf('No file...Quitting!\n');
elseif size(datafile,1) > 1
    checkDIR =uigetfile('.mat');
    load(checkDIR);
    file = checkDIR;
else
    load(datafile(1).name);
    file = datafile(1).name;
end

stim = data.stim;


% SET CALIB FILE HERE
calib = data.FPL.FPLearData;
res.calib = calib;

cd(cwd);

%Call artifact rejection function
[stim] = artifact_rejection(stim);
%Compare GREATtrials to stim.rms-->NEW STIM.RMS
for L = 1:stim.nLevels
    for T = 1:stim.nTrials
        if stim.reject(L,T)
            stim.resp(L,T,:,:) = NaN;
        end
    end
end
%% setting colors
% Colorblind friendly continuous hue/sat changes
cols = [103,0,31;
    178,24,43;
    214,96,77;
    244,165,130;
    253,219,199;
    247, 247, 247;
    209,229,240;
    146,197,222;
    67,147,195;
    33,102,172;
    5,48,97];
cols = cols(end:-1:1, :)/255;

%% Analysis loop


res = MEMRbyLevel(stim);

figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 8 5]}; % xcor, ycor, xwid, yheight

figure;
set(gcf,figure_prop_name,figure_prop_val);
if stim.fc == 7000
    title([subj ' | MEMR - HP | ' condition], 'FontSize', 14)
else
    title([subj ' | MEMR - WB | ' condition], 'FontSize', 14)
end
subplot(1,3,1:2)
semilogx(res.freq / 1e3, res.MEM, 'linew', 2);
xlim([0.2, 8]);
xticks([0.25, 0.5, 1, 2, 4, 8])
%axes('ColorOrder',cols);
%set(gca, 'XTick', ticks, 'XTickLabel', num2str(ticks'), 'FontSize', 14);
legend(num2str(res.elicitor'), 'FontSize', 10, 'NumColumns', 2, 'location', 'best');
xlabel('Frequency (kHz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Ear canal pressure (dB re: Baseline)', 'FontSize', 14, 'FontWeight', 'bold');

power = mean(abs(res.MEM(:, res.ind)), 2);
deltapow = power - min(power);
subplot(1,3,3)

%plot(res.elicitor, mean(abs(res.MEM(:, res.ind)), 2)*5, 'ok-', 'linew', 2);
plot(res.elicitor, deltapow, 'ok-', 'linew', 2);
xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\Delta Absorbed Power (dB)', 'FontSize', 14, 'FontWeight', 'bold');
ymax = max(deltapow+.05);
ylim([0,ymax])
set(gca, 'XScale', 'log', 'FontSize', 14)
drawnow;

res.threshold = interp1(deltapow, res.elicitor, 0.1)

data.res = res; 
%% Export:
cd(datapath);
cd ..
cd('Processed')
if stim.fc == 7000
    fname = [subj,'_MEMR_HP_',condition,  file(end-24:end-4)];
else
    fname = [subj,'_MEMR_WB_',condition,  file(end-24:end-4)];
end

print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'data')
cd(cwd);




