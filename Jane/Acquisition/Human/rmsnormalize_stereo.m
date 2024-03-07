function y = rmsnormalize_stereo(x)
% Returns sound vector rmsnormalized to 0.1
r1 = rms(x(:, 1));
r2 = rms(x(:, 2));
y = x*0.02/max(r1, r2);