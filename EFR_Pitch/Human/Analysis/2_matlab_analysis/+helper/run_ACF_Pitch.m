function [win_ACFs,t_acf,P1,P2] = run_ACF_Pitch(efrs,F0,fs,plot_acg)
%RUN_ACF_PITCH - input = efr matrix (efrs), Fundamental Freq (F0)
% and sample rate (fs)
% sorted by harmonic rank 
%output is the windowed ACFs, and time in ms, P1 (periodicity at  and P2
%TODO: implement this method in the main scripts...save space lol

if nargin == 3
    plot_acg = 0;
end

if plot_acg
    acf_fig = figure();
    load('cmap.mat');
end

window = 1.5/F0; %window in seconds
%Filter EFRs
[b,a] = butter(6,[60,500]/(fs/2));
efrs = filtfilt(b,a,efrs);
acf_bin_length = round(window*fs);
inds = 1:acf_bin_length/10:(length(efrs)-round(window*fs));
inds = round(inds);

%I bet there's a more efficient way to do this...
for i = 1:size(efrs,2)
    for j = 1:length(inds-1)
        chunk = efrs(inds(j):(inds(j)+acf_bin_length-1),i); %temp
        chunk_corr = xcorr(chunk,'normalized');
        temp = chunk_corr(acf_bin_length:end);
        acf_out(j,:) = temp.*(temp>0);
    end
    
    if plot_acg
        subplot(1,size(efrs,2),i)
        surf((inds)/fs, 1e3*(1:size(acf_out,2))/fs, acf_out','EdgeColor','interp');
        colormap(cmap)
        view(2);
        title(['Response ',num2str(i)]);
    end

    win_ACFs(i,:) = mean(acf_out);
end

pk1_t = (1/F0)*fs;
pk2_t = (1/(2*F0))*fs;

P1 = win_ACFs(:,round(pk1_t));
P2 = win_ACFs(:,round(pk2_t));

t_acf = 1000*(1:size(acf_out,2))/fs;

%plotting
if plot_acg
    han=axes(acf_fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Time-lag (ms)');
    xlabel(han,'Time(s)');
    acgtitle = strcat('Autocorrelogram');
    sgtitle(acgtitle)
end

end

