if ~exist('presentedSNRs', 'var')
    error('Please load a file and then run this code');
end
[up, down, all] = getReversals(presentedSNRs);
threshold = mean(presentedSNRs(up((end-3):end)))*0.25 + ...
    0.75*mean(presentedSNRs(down((end-3):end)))

plot(presentedSNRs, 'o-k', 'linew', 2, 'MarkerSize', 10);
hold on;
plot(1:numel(presentedSNRs), ones(size(presentedSNRs))*threshold, ...
    'r--', 'linew', 2);
ylabel('SNR (dB)', 'FontSize', 14);
xlabel('Trial Number', 'FontSize', 14);

