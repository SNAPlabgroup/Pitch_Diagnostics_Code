%% Analysis Script for analyzing MRT data 
% Author: Samantha Hauser with code from Hari Bharadwaj
% Created: August 2023
% Last Updated: August 17, 2023
% Purpose:
% Helpful info:

%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,[ subj, '_block13*.mat'])).name};

numOfFiles = length(datafile);
if numOfFiles > 1
    fprintf('More than one file...Select the right one?\n');
    datafile = uigetfile([ subj, '_block13*.mat']); 
elseif numOfFiles < 1
    fprintf('No files for this subject...Quitting.\n')
    cd(cwd); 
    return
else
    datafile = datafile{1}; 
end 

load(datafile);

cd(cwd);

%% Colors
red = [[254,229,217]/255; [252,187,161]/255; [252,146,114]/255; [251,106,74]/255; [222,45,38]/255; [165,15,21]/255]; 

%% From plotResults.m 

SNRs = -20:4:4; % OLD: -20:5:5;
ntrials = zeros(size(SNRs)); 
score = zeros(size(SNRs)); 
scorestd = zeros(size(SNRs)); 

for k = 1:numel(SNRs)
    resp_SNR = responseTable(responseTable(:, 4) == SNRs(k), 2);
    target_SNR = responseTable(responseTable(:, 4) == SNRs(k), 3);
    ntrials(k) = sum(responseTable(:, 4) == SNRs(k)); 
    score(k) = sum(resp_SNR == target_SNR) / ntrials(k); 
    scorestd(k) = score(k)*(1-score(k))/sqrt(ntrials(k)); 
end

%% Sigmoid fit
ft = fittype('L + (U-L)/(1 + exp(-4*log(3)*(x-xmid)/xscale80))','indep','x') 
mdl = fit(SNRs',score',ft, 'start', [0.4, .9, 0,1]) 
fitx = SNRs(1)-4:.1:SNRs(end)+4; 
fity = mdl(fitx); 
ind = find(fity > .7599,1,  'first'); 
thresh = fitx(ind); 
%% Plotting
figure; 
errorbar(SNRs - 0.5, score, scorestd, 'sr', 'linew', 2);
hold on; 
plot(fitx, fity, 'k', 'linew', 2); 
plot([-24, thresh], [.76, .76], 'o--k', 'linew', 2)
plot([thresh, thresh], [.2, .76], 'o--k', 'linew', 2)
xlabel('SNR (dB)', 'fontsize', 14, 'FontWeight', 'bold');
ylabel('Proportion Correct', 'fontsize', 14, 'FontWeight', 'bold');
text(thresh+1, .25, sprintf('%.2f dB', thresh), 'FontSize', 14); 
xlim([-24, 4]); 
ylim([.2, 1]); 
grid on;
set(gca, 'FontSize', 14, 'XTick', SNRs);
title([subj ' | MRT Tone Complex | ' condition], 'FontSize', 14)
set(gca, 'FontSize', 14)
drawnow; 

res.ntrials = ntrials; 
res.score = score; 
res.scorestd = scorestd;
res.SNRs = SNRs; 
res.fitx = fitx; 
res.fity = fity; 
res.thesh76 = thresh; 
  
%% Save Variables and figure
cd(datapath);
fname = [subj,'_MRT_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'res')
cd(cwd);