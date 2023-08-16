function plotTuningCurve(subj)

if nargin < 1
    subj = 'S363'; % Enter subj to look at here
end

folder = ['C:\Experiments\Sam\TuningCurveTDT\Results\', subj];
files = dir([folder, '\', subj, '_Tuning_*.mat']);

if length(files) < 3
    fprintf('Only %d trials. Run %d more.\n', length(files), 3-length(files));
    done = false;
end

for i = 1:length(files)
    res = load([files(i).folder, '\', files(i).name]);
    
    % [FileName,PathName,FilterIndex] = uigetfile(strcat('./Results/*_Tuning_*.mat'),...
    %     'Please pick MEM data file to analyze');
    % TCfile = fullfile(PathName, FileName);
    % res = load(TCfile);
    
    figure(92);
    
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
    semilogx(fadj/1000,Lplot,'k','linewidth',1);
    set(gca, 'XTick', [2, 2.5, 3, 3.5, 4, 5]);
    grid on;
    xlabel('Frequency (kHz)');
    ylabel('Threshold (dB SPL)');
    set(gca,'fontsize',16);
    
    f(:,i) = fadj;
    L(:,i) = Lplot;
end

if done
    figure(92);
    hold on;
    semilogx(mean(f/1000, 2),mean(L,2),'r','linewidth',2);
end
end