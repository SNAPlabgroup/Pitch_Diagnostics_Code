[FileName,PathName,FilterIndex] = uigetfile(strcat('*_Tuning_*.mat'),...
    'Please pick MEM data file to analyze');
TCfile = fullfile(PathName, FileName);
res = load(TCfile);

ftrack = res.results.ftrack;
indftrack = 700:5000; % Hardcoding for roughly 2kHz to 5 kHz

ftrack  = ftrack(indftrack);
Ltrack = res.results.Ltrack(indftrack);

% Register based on the tip and adjust for response delay shift
[Lmin, imin] = min(Ltrack);
fmin = ftrack(imin);
fadj = ftrack * 4000/fmin;

% Smoothing
framelen = 999;
wts = dpss(framelen, 1, 1);
order = 3;
Lplot = sgolayfilt(Ltrack, order, framelen, wts);

hold on;
semilogx(fadj/1000,Lplot,'k','linewidth',2);
set(gca, 'XTick', [2, 2.5, 3, 3.5, 4, 5]);
grid on;
xlabel('Frequency (kHz)');
ylabel('Threshold (dB SPL)');
set(gca,'fontsize',16);

