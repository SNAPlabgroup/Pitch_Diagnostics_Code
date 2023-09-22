function stim = makeCABRstim(clickdur, rate, dur, fs, jitterwindow, polarity, plotclickspec)
%  Make condensation click pulse sequences.
%
% USAGE:
%  stim = makeCABRstim(clickdur, rate, dur, fs, jitterwindow, polarity, plotclickspec)
%
% Example: stim = makeCABRstim(clickdur, 39, 3, 48828.125, 0.005);
%
% clickdur: duration of click (micro.s): will be rounded to integer samples
% rate: Average click rate (also see jitterwindow) (Hz)
% dur: Total stimulus duration (s)
% fs: Sampling rate (Hz)
% jitterwindow: The width of the time range over which to jitter ICI.
%               (uniformly). Note that the center of the ICI distribution
%               is 1/rate. For example, if rate is 10, and jitterwindow is
%               0.05, then the ICI is distributed uniformly over [75,
%               125] ms.
% polarity: 1 or -1 for condensation or rarefaction clicks
% plotclickspec: If non-zero, plot the spectrum of a single click (default 0)
%----------------------------------------------------------------------
% Hari Bharadwaj, Apr 25, 2015
% hari@nmr.mgh.harvard.edu
%----------------------------------------------------------------------

if(~exist('plotclickspec', 'var'))
    plotclickspec = 0;
end
single_click = zeros(1, floor((1/rate - jitterwindow/2)*fs));
nsamps_single = numel(single_click);
clickdur_samps =  floor(clickdur*1e-6*fs);
clickdur_ind = 0:(clickdur_samps-1);
clickdur_ind = clickdur_ind - round(median(clickdur_ind));
if polarity~=0
    polarity = sign(polarity);
else
    error('Polarity cannot be zero!');
end
single_click(floor(nsamps_single/2) + clickdur_ind) = 0.95*polarity;
nclicks = ceil(dur*rate)+ 50; % Just some cushion to index out 'dur' seconds
jitters = round(rand(1, nclicks) * jitterwindow * fs);
stim = single_click;
for k = 1:nclicks
    stim = [stim, zeros(1, jitters(k)), single_click]; %#ok<AGROW>
end

stim = stim(1:floor(dur*fs)); % Here is where the cushion comes handy

if(plotclickspec ~= 0)
    [pxx, f] = pmtm(single_click, 2, [], fs);
    figure;
    plot(f*1e-3, pow2db(pxx), 'linew', 2);
    xlabel('Frequency (kHz)', 'FontSize', 16);
    ylabel('Power (dB)', 'FontSize', 16);
end







