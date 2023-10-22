fs = 48828; % Sampling rate
duration = 1.0; % Second
fm = 223; % Modulation frequency
fc = 4000; % Carrier frequency

t = 0:(1/fs):(duration - 1/fs);
carr = sin(2*pi*fc*t);

mod = (square(2*pi*fm*t, 25) + 1)/2;
signal = mod .* carr;

audiowrite('RAM4k.wav',signal,fs)