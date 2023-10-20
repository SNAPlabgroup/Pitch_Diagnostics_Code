clear x;
nruns = numel(allResults);
x = 0;
f = 1500:5:5600;
for index = 1:nruns
    % Which run of the loaded data to plot?
    temp = allResults{index};
    x = x + interp1(temp(:, 1), temp(:, 2), f);
    x_std(index,:)= interp1(temp(:, 1), temp(:, 2), f); 
end
x = x/nruns;
hold on;
plot(f, x, 'k');
hold on;
nfilt = 99; % How much to smooth
ord = 3; % filter order
plot(f,...
    flipud(sgolayfilt(flipud(sgolayfilt(x, ord, nfilt)), ord, nfilt)),...
    'r', 'linew', 2);
set(gca,'XLim',[1600, 5400], 'fontsize', 16,...
    'xscale', 'log', 'box', 'off');
ylabel('Masking Threshold (dB SPL)', 'fontsize', 16);
xlabel('Noise Center Frequency (Hz)', 'fontsize', 16);
legend('Raw data', 'Smoothed data');
grid on;

