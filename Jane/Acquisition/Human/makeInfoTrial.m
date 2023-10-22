function x = makeInfoTrial(T, M1, M2, inds, loc, TMR, staggerT, staggerM, fs)
% Make the stimulus for a single trial of the informational masking
% speech-on-speech paradigm.
%
% USAGE:
%   [x, fs] = makeInfoTrial(T, M1, M2, inds, loc, fs);
%
% INPUTS:
%   T - Target voice (string: one from 1-12M and 1-12F. e.g., '12F')
%   M1 - First masker voice (presented on left, if loc ~= 0)
%   M2 - Second masker voice (presented on right, if loc ~= 0)
%   inds - matrix of 3 rows and 5 columns saying what the index of words
%       is from the reference cell array. First row is for target, second
%       is for M1 and thrid is for M2.
%   loc - Either 0, 15 or 45 degree separation
%   TMR - TMR in dB (default is 0)
%   staggerT - Stagger after target onset (seconds, default is 0.3);
%   staggerM - Inter-masker onset stagger (seconds, default is 0.3);
%   fs - desired sampling rate (Hz, default is 48828.125)
%
% OUTPUTS:
%   x - The stereo sound (length varies depending on choice of talkers and
%       words.
%
% REFERENCES:
%   Kidd Jr, G., Best, V., & Mason, C. R. (2008). Listening to every other
%   word: Examining the strength of linkage variables in forming streams
%   of speech. J Acoust Soc Am, 124(6), 3793-3802.
%
%--------------------------------------------------------------------------
% Copyright 2019 Hari Bharadwaj. All rights reserved.
% hbharadwaj@purdue.edu
% -------------------------------------------------------------------------

if ~exist('TMR', 'var')
    TMR = 0;
end

if ~exist('staggerT', 'var')
    staggerT = 0.3; 
end

if ~exist('staggerM', 'var')
    staggerM = 0.3; % Inter-masker onset stagger
end

if ~exist('fs', 'var')
    fs = 48828.125;
end

% Load word matrix -- provides the variable 'words'
load('C:\Experiments\SIN_INFO\BUC-Mx_2014_FiveCategories\BUC_Mx_2014_Text.mat');

audiodir = 'C:\Experiments\SIN_INFO\BUC-Mx_2014_FiveCategories\';
% Get target audio
for k = 1:5
    % Target is always at the center
    fname = strcat(audiodir, '0_BU\0_BU_', words{inds(1, k), k}, '_', T, '.wav'); %#ok<USENS>
    [temp, fs_wav] = audioread(fname);
    [P, Q] = findresamplePQ(fs, fs_wav);
    Taudio{k} = rmsnormalize_stereo(resample(temp, P, Q));
end

% Get M1 audio (on the left if loc ~= 0)
for k = 1:5
    switch loc
        case 0
            fname = strcat(audiodir, '0_BU\0_BU_', words{inds(2, k), k}, '_', M1, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M1audio{k} = rmsnormalize_stereo(resample(temp, P, Q)); %#ok<*AGROW>
        case 15
            fname = strcat(audiodir, 'L15_BU\L15_BU_', words{inds(2, k), k}, '_', M1, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M1audio{k} = rmsnormalize_stereo(resample(temp, P, Q));
        case 45
            fname = strcat(audiodir, 'L45_BU\L45_BU_', words{inds(2, k), k}, '_', M1, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M1audio{k} = rmsnormalize_stereo(resample(temp, P, Q));
            
    end
end

% Get M2 audio (on the right if loc ~= 0)
for k = 1:5
    switch loc
        case 0
            fname = strcat(audiodir, '0_BU\0_BU_', words{inds(3, k), k}, '_', M2, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M2audio{k} = rmsnormalize_stereo(resample(temp, P, Q));
        case 15
            fname = strcat(audiodir, 'R15_BU\R15_BU_', words{inds(3, k), k}, '_', M2, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M2audio{k} = rmsnormalize_stereo(resample(temp, P, Q));
        case 45
            fname = strcat(audiodir, 'R45_BU\R45_BU_', words{inds(3, k), k}, '_', M2, '.wav');
            [temp, fs_wav] = audioread(fname);
            [P, Q] = findresamplePQ(fs, fs_wav);
            M2audio{k} = rmsnormalize_stereo(resample(temp, P, Q));
            
    end
end

%% Stitch together words to form stimulus



for k = 1:5
    len1 = length(M1audio{k});
    len2 = length(M2audio{k});
    if len1 > len2
        pad = zeros(len1 - len2, 2);
        M2audio{k} = [M2audio{k}; pad];
    else
        pad = zeros(len2 - len1, 2);
        M1audio{k} = [M1audio{k}; pad];
    end
    
    leader = 1 + round(rand); % Should be 1 or 2
    pad = zeros(ceil(staggerM*fs), 2);
    if leader == 1
        M2audio{k} = [pad; M2audio{k}];
        M1audio{k} = [M1audio{k}; pad];
    else
        M2audio{k} = [M2audio{k}; pad];
        M1audio{k} = [pad; M1audio{k}];
    end
    % Combined masker audio
    Maudio{k} = M1audio{k} + M2audio{k};
end

% Concatenate target and maskers at the right TMR
factor = db2mag(TMR);
GAP = ceil((2*staggerT + staggerM)*fs);
Tons = (0:4)*GAP;
Mons = Tons + ceil(staggerT*fs);

x = zeros(GAP*6, 2);
for k = 1:5
    x(Tons(k) + (1:length(Taudio{k})), :) = x(Tons(k) + (1:length(Taudio{k})), :) + Taudio{k}*factor;
    x(Mons(k) + (1:length(Maudio{k})), :) = x(Mons(k) + (1:length(Maudio{k})), :) + Maudio{k};
end

