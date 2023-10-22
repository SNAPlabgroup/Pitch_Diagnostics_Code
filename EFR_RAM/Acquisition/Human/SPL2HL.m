function HL = SPL2HL(SPL, fc)
% Converts dB SPL to HL assuming SNAPlab Sennheiser HDA 300 calibrations


if any(fc < 125) || any(fc > 16000)
    error('Only frequencies between 125 Hz and 16 kHz are valid!');
end
f = [125, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000,...
    8000, 9000, 10000, 11200, 12500, 14000, 16000];
        
offset = [26.2, 20.1, 8.6, 5.1, 2.7, 3.2, 0.5, -1.6, 0.1, 11.3,...
    20.9, 23.1, 27.1, 18.5, 22.9, 27.0, 32.8, 47.7];

off = interp1(f, offset, fc, 'spline');

HL = SPL - off;