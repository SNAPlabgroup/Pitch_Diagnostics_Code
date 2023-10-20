N = [31, 50, 50, 50, 50, 50, 31]; % Number of trials per SNR (variable)
SNRs = 4:-4:-20; % dB

% Setting random generator seed and state
load('s.mat');
rng(s);

nconds = numel(SNRs);
nwordsperlist = 6;
totaltrials = sum(N);
targets = repmat(1:nwordsperlist, 1, (totaltrials/nwordsperlist));

targets = targets(randperm(totaltrials));

% Exclude M5 because of different naming convention of files and to balance
% male and female voices.
speakers = {'F1', 'F2', 'F3', 'F4', 'M1', 'M2', 'M3', 'M4'};
speakers = speakers(randi(numel(speakers), [1, totaltrials]));
wordlists = 1:50;
wordlists = wordlists(randi(numel(wordlists), [1, totaltrials]));

SNRcounts = zeros(nconds, 1);
SNRlist = [];
for kk = 1:max(N)
    for k = 1:nconds
        if SNRcounts(k) < N(k)
            SNRlist = [SNRlist, SNRs(k)]; %#ok<AGROW>
            SNRcounts(k) = SNRcounts(k) + 1;
        end
    end
end
for k = 1:totaltrials
    target = targets(k);
    wordlist = wordlists(k);
    fname = [speakers{k}, '_b', num2str(wordlist), '_w',...
        num2str(target), '_orig.wav'];
    fpath = ['./audio/', speakers{k}];
    fnamefull = fullfile(fpath, fname);
    
    % Load file and process
    x = resample(audioread(fnamefull), 4069, 4000);
    fs = 48828;
    SNR = SNRlist(k);
    y = stonemoore2014(x, fs, SNR);
    
    % Save processed audio
    savename = ['./trialaudio/trial', num2str(k), '.mat'];
    save(savename, 'y', 'fs', 'SNR', 'target', 'wordlist');
    if mod(k, 24) == 0
        fprintf(1, '###########\n Done with %d / %d trials ...\n#########\n', k, totaltrials);
    end
end
    
