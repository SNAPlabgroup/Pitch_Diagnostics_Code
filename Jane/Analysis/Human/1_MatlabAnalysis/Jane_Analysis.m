%% Analysis Script for analyzing Jane Task Data
% Author: Samantha Hauser with code from Hari Bharadwaj
% Created: August 2023
% Last Updated: August 17, 2023
% Purpose:
% Helpful info:

%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,[ subj, '_*.mat'])).name};
cd(cwd);
numOfFiles = length(datafile);
if numOfFiles == 0
    fprintf('No Files for this subject\n')
    return 
end
threshold = zeros(numOfFiles, 1); 

%% Colors
red = [[254,229,217]/255; [252,187,161]/255; [252,146,114]/255; [251,106,74]/255; [222,45,38]/255; [165,15,21]/255]; 

%% Loop for each block
for i = 1:numOfFiles
    
    load([datapath, filesep, datafile{i}]);
    
    [up, down, all] = getReversals(presentedSNRs);
    threshold(i,1) = mean(presentedSNRs(up((end-3):end)))*0.25 + ...
        0.75*mean(presentedSNRs(down((end-3):end))); 
    
    plot(presentedSNRs, 'o-k', 'linew', 2, 'MarkerSize', 10 ); % 'Color', [150, 150, 150]./255
    hold on;
    plot(1:numel(presentedSNRs), ones(size(presentedSNRs)).*threshold(i,1), ...
        '--', 'linew', 1.5, 'color', red(i, :));

end
res.thresholds = threshold; 
res.threshold_avg = mean(threshold); 
res.threshold_std = std(threshold); 
res.numOfBlocks = numOfFiles; 

hold on; 
plot(1:45, ones(1, 45).*res.threshold_avg, ...
        '-', 'linew', 2, 'color', red(6,:));
plot(1:45, ones(1, 45).*(res.threshold_avg + res.threshold_std), ...
        ':', 'linew', 2, 'color', red(6,:));    
plot(1:45, ones(1, 45).*(res.threshold_avg - res.threshold_std), ...
        ':', 'linew', 2, 'color', red(6,:));  
text(35, res.threshold_avg-res.threshold_std-4, sprintf('%.2f dB', res.threshold_avg), 'FontSize', 14);  
    
title([subj ' | Jane Task | ' condition], 'FontSize', 14)
set(gca, 'FontSize', 14)
%xlim([0, 50])
%ylim([-50, 15])
ylabel('SNR (dB)', 'FontWeight', 'bold')
xlabel('Trial Number', 'FontWeight', 'bold')
legend('Trial', 'Threshold')
drawnow; 
hold off; 

%% Save Variables and figure
cd(datapath);
fname = [subj,'_Jane_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'res')
cd(cwd);