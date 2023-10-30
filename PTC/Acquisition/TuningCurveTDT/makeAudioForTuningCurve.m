[noise, fs, startF, endF, bw, freqseq] = makeNoiseSweep(2);

pulse = maketone(0, 0, 4000, 0.5, 0.02, fs, 0.2)';

Npulses_before_noise = 5;

pulsedur_samps = numel(pulse);

noise = [zeros(1, pulsedur_samps * Npulses_before_noise), noise];
freqseq = [zeros(1, pulsedur_samps * Npulses_before_noise), freqseq];
nreps = ceil(numel(noise)/pulsedur_samps);

tone = [];
for k = 1:nreps
    tone = [tone, pulse]; %#ok<AGROW> 
end

excess_samps = numel(tone) - numel(noise);
noise = [noise, zeros(1, excess_samps)];
freqseq = [freqseq, zeros(1, excess_samps)];
save('tuningStims', 'fs', 'freqseq', 'noise', 'tone', 'bw', 'startF', 'endF');
