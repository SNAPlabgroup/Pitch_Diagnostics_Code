function [P, Q] = findresamplePQ(f_desired, f_orig)
% Finds the upsampling and downsampling factors to get approximate desired
% sampling rates. Sort of finding rational approximation.
% USAGE:
%   [P, Q] = findresamplePQ(f_desired, f_orig);

a = factor(round(f_desired));
b = factor(round(f_orig));

commonfac = prod(a(ismember(a, b)));

P = round(f_desired) / commonfac;
Q = round(f_orig) / commonfac;