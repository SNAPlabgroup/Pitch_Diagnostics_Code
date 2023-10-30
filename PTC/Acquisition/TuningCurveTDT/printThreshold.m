%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN AFTER LOADING SOME RESULT FILE MANUALLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramlist = responseList(:, 4);
revList = [];
downList = [];
upList = [];
nReversals = 0;
for k = 4:numel(paramlist)
    if((paramlist(k-1) > paramlist(k)) && (paramlist(k-1) == paramlist(k-2)) ...
            && (paramlist(k-2) > paramlist(k-3)))
        nReversals = nReversals + 1;
        revList = [revList, (k-1)]; %#ok<*AGROW> 
        downList = [downList, (k-1)];
    end
end

for k = 3:numel(paramlist)
    if((paramlist(k-1) < paramlist(k)) && (paramlist(k-1) < paramlist(k-2)))
        nReversals = nReversals + 1;
        revList = [revList, (k-1)];
        upList = [upList, (k-1)];
    end
end

for k = 4:numel(paramlist)
    if((paramlist(k-1) < paramlist(k)) ...
            && (paramlist(k-1) == paramlist(k-2)) ...
            && (paramlist(k-2) < paramlist(k-3)))
        nReversals = nReversals + 1;
        revList = [revList, (k-1)];
        upList = [upList, (k-1)];
    end
end

revList = sort(revList, 'ascend');
select = paramlist(revList(end-6:end));

%% Should be the same in trial audio and here for accurate thresholds
minparam = 1;
maxparam = 40;
params = minparam:maxparam;
Ls = linspace(-20, 100, numel(params));

thresh_mean = mean(Ls(select));
thresh_std = std(Ls(select))/sqrt(7);

fprintf(['Threshold = %0.0f ', char(177),...
    ' %0.1f dB\n'],...
    thresh_mean, thresh_std);

figure;
plot(Ls(paramlist), 'k-', 'linew', 2);
hold on;
yline(thresh_mean, 'r-','LineWidth', 2);
yline(thresh_mean + 1.96 * thresh_std, 'r--', 'LineWidth', 2);
yline(thresh_mean - 1.96 * thresh_std, 'r--', 'LineWidth', 2);

set(gca, 'FontSize', 16);
xlabel('Trial #', 'FontSize', 16);
ylabel(['Tone level (', char(956), 'dB re: max)'], 'FontSize',16);
legend('Adaptive Track', 'Estimated Threshold', ...
   '95% CI');
text(numel(paramlist)/2, mean(Ls(paramlist)),...
    [num2str(thresh_mean), ' dB'] ,'FontSize', 16);