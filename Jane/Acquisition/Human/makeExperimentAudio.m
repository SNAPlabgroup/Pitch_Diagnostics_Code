SNRlist = 18:-3:-54;
Nper = 20;
T = '2F';
M1 = '1F';
M2 = '6F';
loc = 15;
stagger = 0.325;
fs = 48828;
audiodir = 'C:\Experiments\SIN_INFO\trialaudio\';
for kSNR = 1:numel(SNRlist)
    SNR = SNRlist(kSNR);
    fprintf(1, 'Creating files for SNR = %d ...\n', SNR);
    for k = 1:Nper
        fname = strcat(audiodir, 'INFOMASK_SNR', num2str(SNR),...
            '_trial_', num2str(k), '.mat');
        inds = makeinds;
        x = makeInfoTrial(T, M1, M2, inds, loc, SNR, stagger, stagger, fs);
        save(fname, 'x', 'inds');
    end
end

        

