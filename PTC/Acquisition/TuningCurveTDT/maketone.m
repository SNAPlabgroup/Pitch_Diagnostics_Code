function x = maketone(L, Lref, fc, dur, rise, fs, pad)
% Make a tone of a specified level (relative to a fullscale reference), 
% frequency, duration, rise time and post-padding.
%
% USAGE:
%   x = maketone(L, Lref, fc, dur, rise, fs, pad);
% 
% INPUTS:
%   L - Tone level
%   Lref - Level expected for a full-scale (-1 to 1) tone
%   fc - Center frequency of tone (Hz)
%   dur - Duration of tone (seconds)
%   rise - Rise time of tone (seconds)
%   fs - Sampling rate
%   pad - Duration of post padding of zeros. Makes it easier to construct
%   sequences
%
% OUTPUT:
%   x - Tone with the desired parameters (column vector)
% ---------------------------------------------------
% Copyright 2023 Hari Bharadwaj. All rights reserved.
% hari.bharadwaj@pitt.edu
% ---------------------------------------------------

t = 0:(1/fs):(dur - 1/fs);

A = db2mag(L - Lref);

tone = A * rampsound(sin(2 * pi * fc * t(:)), fs, rise);
x = [tone; zeros(ceil(fs * pad), 1)];


