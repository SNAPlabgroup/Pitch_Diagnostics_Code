%% Use this general template to generate stimuli for the adaptive task
% Regardless of the actual units/values, the stimulus parameter is varied
% over a grid of from minparam (set to 1), and maxparam (set to 40).

Nfiles = 5; % Number of files per parameter level
minparam = 1;
maxparam = 40;
params = minparam:maxparam;

% Here, the parameter being adapted is tone level in dB units
% If the typical unit is linear, it would generally be appropriate to use
% 'logspace' here. 
Ls = linspace(-20, 100, numel(params));

% Other stimulus parameters are set here
dur = 250e-3;
pad = 400e-3;
fs = 48828.125;
rise = 20e-3;
Lref = 114;  % ER-2 dB SPL for a full-scale (i.e., RMS 1/sqrt(2)) tone.
fc = 4000;
Ldummy = -Inf;

for k = 1:numel(params)
    for nf = 1:Nfiles
        answer = randi(3);

        L = Ls(k);
        param = params(k);

        dummy = maketone(Ldummy, Lref, fc, dur, rise, fs, pad);
        target = maketone(L, Lref, fc, dur, rise, fs, pad);
        switch answer
            case 1
                y = [target; dummy; dummy];
            case 2
                y = [dummy; target; dummy];
            case 3
                y = [dummy; dummy; target];
        end
        fname = strcat('./trialaudio/trial', num2str(param),...
            '_', num2str(nf), '.mat');
        save(fname, 'y', 'fs', 'param', 'answer');
    
    end
end

