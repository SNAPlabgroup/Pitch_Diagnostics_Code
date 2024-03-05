% SFOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: October 27, 2023
% Purpose:
% Helpful info: Need to add Qerb calculation and consider eliminating long
% components (w/ IFFT method)

%%% Set these parameters %%%%%%%%%%%%%

windowdur = 0.038; % 40ms in paper
offsetwin = 0.0; % 20ms in paper
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
cwd = pwd;
cd([datapath, filesep, 'Preprocessed'])
datafile = dir(fullfile(cd,['sweptSFOAE_*.mat']));
if length(datafile) < 1
    fprintf('No file...Quitting!\n');
elseif size(datafile,1) > 1
    checkDIR =uigetfile('.mat');
    load(checkDIR);
    file = checkDIR; 
else
    load(datafile(1).name);
    file = datafile(1).name; 
end

cd(cwd);

%% Get data structure and handle old data storage formats
stim = data.stim;

if ~isfield(stim, 'scale')
    stim.scale = 'log';
    stim.nearfreqs = [1.1,1.12, 1.14,1.16];
end

figure; plot(stim.SuppBuffs(1,1:400)); hold on; plot([128, 247], [0,0], 'or')
text(128, .1, '128'); text(247, .1, '247')
ask_delay = inputdlg('extra delay?');  % 247; %128
delay_oops = str2double(ask_delay{1}); 


% Check Awake or Sedated
sed = questdlg(sprintf('Awake or Sedated? Dir name is: %s', data.info.dir), 'Sedated?', 'Awake', 'Sedated', 'Awake'); 
switch sed
    case 'Sedated'
        sedated_flag = 1; 
    case 'Awake'
        sedated_flag = 0; 
end
%% Get appropriate calibration file
calib = data.FPL.FPLearData;
res.calib = calib; 

%% Analysis Set up
% Set variables needed from stim.
phiProbe_inst = 2*pi*stim.phiProbe_inst;
t = stim.t;

% downward vs upward sweeps
if stim.speed < 0
    f1 = stim.fmax;
    f2 = stim.fmin;
else
    f1 = stim.fmin;
    f2 = stim.fmax;
end

% set freq we're testing and the timepoints when they happen.
if strcmp(stim.scale, 'log')    %linear sweep
    testfreq = 2 .^ linspace(log2(f1), log2(f2), npoints);
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else
    testfreq = linspace(f1, f2, npoints);
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end

%duration changes w/ frequency
%durs = .038*(2.^(-0.3*(t_freq-stim.buffdur))-1)/ (-0.3*log(2)) + 0.038;
durs = 0.0244.*(t_freq-stim.buffdur)+0.038; %linear change

nfreqs = stim.nearfreqs;
%% Artifact rejection

% Cancel out stimulus
SFOAEtrials = zeros(size(stim.ProbeBuffs)); 
subtraction = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs;
SFOAEtrials(:,1:end-delay_oops) = subtraction(:,delay_oops +1:end);

trials = size(SFOAEtrials,1);

% high pass filter the response (can also be done on ER10X hardware)
filtcutoff = 300;
b = fir1(1000, filtcutoff*2/stim.Fs, 'high');
SFOAEtrials= filtfilt(b, 1, SFOAEtrials')';

% Set empty matricies for next steps
coeffs = zeros(npoints, 2);
a_temp = zeros(trials, npoints);
b_temp = zeros(trials, npoints);

% Least Squares fit of SF Only for AR
for x = 1:trials
    SFOAE = SFOAEtrials(x, :);
    fprintf(1, 'Checking trial %d / %d for artifact\n', x, trials);
    
    for k = 1:npoints
        windowdur = durs(k);
        win = find( (t > (t_freq(k) - windowdur/2)) & ...
            (t < (t_freq(k) + windowdur/2)));
        taper = hanning(numel(win))';
        model_sf = [cos(phiProbe_inst(win)) .* taper;
            -sin(phiProbe_inst(win)) .* taper];
        resp = SFOAE(win) .* taper;
        coeffs(k, 1:2) = model_sf' \ resp';
    end
    a_temp(x,:) = coeffs(:, 1);
    b_temp(x,:) = coeffs(:, 2);
end

oae = abs(complex(a_temp, b_temp));
median_oae = median(oae);
std_oae = std(oae);

resp_AR = SFOAEtrials;
for j = 1:trials
    for k = 1:npoints
        if oae(j,k) > median_oae(1,k) + 3*std_oae(1,k) %3*SD to match DP script
            win = find( (t > (t_freq(k) - durs(k).*.1)) & ...
                (t < (t_freq(k) + durs(k).*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end

SFOAE = mean(resp_AR, 'omitNaN'); % mean SFOAE after artifact rejection

%% LSF Analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
tau = zeros(npoints, 1);
coeffs_noise = zeros(npoints,8);

% Generate model of chirp and test against response
for k = 1:npoints
    fprintf(1, 'Running window %d / %d\n', k, (npoints));
    windowdur = durs(k);
    win = find( (t > (t_freq(k)-windowdur/2)) & ...
        (t < (t_freq(k)+windowdur/2)));
    taper = hanning(numel(win))';
    
    % set the models
    
    % SF probe frequency
    model = [cos(phiProbe_inst(win)) .* taper;
        -sin(phiProbe_inst(win)) .* taper];
    
    % nearby frequencies for nf calculation
    model_noise = [cos(nfreqs(1)*phiProbe_inst(win)) .* taper;
        -sin(nfreqs(1)*phiProbe_inst(win)) .* taper;
        cos(nfreqs(2)*phiProbe_inst(win)) .* taper;
        -sin(nfreqs(2)*phiProbe_inst(win)) .* taper;
        cos(nfreqs(3)*phiProbe_inst(win)) .* taper;
        -sin(nfreqs(3)*phiProbe_inst(win)) .* taper;
        cos(nfreqs(4)*phiProbe_inst(win)) .* taper;
        -sin(nfreqs(4)*phiProbe_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 2);
    resid = zeros(maxoffset, 1);
    coeff_noise = zeros(maxoffset, 8);
    
    for offset = 0:maxoffset
        resp = SFOAE(win+offset) .* taper;
        
        coeff(offset + 1, :) = model' \ resp';
        coeff_noise(offset +1, :) = model_noise' \ resp';
        
        resid(offset +1) = sum( (resp - coeff(offset+1, :) * model).^2);
    end
    
    [~, ind] = min(resid);
    
    coeffs(k, :) = coeff(ind, :);
    coeffs_noise(k,:) = coeff_noise(ind,:);
    
    tau(k) = (ind - 1) * (1/stim.Fs); % delay in sec
    
end

%% Amplitude and delay calculations
a = coeffs(:, 1);
b = coeffs(:, 2);

phi = tau.*testfreq'; % cycles (from delay/offset)
phasor = exp(-1j * phi* 2 * pi);

% for noise
noise2 = zeros(npoints,4);
for i = 1:2:8
    noise2(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end

oae_complex = complex(a, b).*phasor;
noise_complex = mean(noise2,2);
res.multiplier = stim.VoltageToPascal.* stim.PascalToLinearSPL;

%% Windowing
% Can use IFFT method from Hari's additon to early DPOAE code
% use energy between 0.5 and 1.5 ms



%% Separating D and R components by IFFT
% Code from IFFT method of separating D and R components from HB. 

% win_center = stim.fmin; 
% win_ind = 1;
% bw = []; 
% 
% while win_center(win_ind) <= stim.fmax
%     % based on Abdala 2022, "Hann window varied as sqrt of frequency from
%     % 400 Hz at 500 to 1600 Hz at 8000"
%     bw(win_ind) = sqrt(win_center(win_ind))*40/sqrt(5);
%     
%     % Hann windowed segments with 50 Hz overlap 
%     win_center(win_ind+1) = win_center(win_ind)+(bw(win_ind)/2)-50;
%     
%     win_ind = win_ind+1; 
% end
% 
% % window each segement of complex oae and nf
% 
% fs = stim.Fs;
% timedur = 30e-3;
% % Check if 1/(frequency resolution) is long enough
% if 1/mean(abs(diff(testfreq))) < timedur/2
%     warning('Too few DPOAE points, aliasing  will likely  occur');
% end
% N = roundodd(timedur * fs);
% f = (0:(N-1))*fs/N; % FFT bin frequencies
% FFT_sf =  interp1(testfreq, oae_complex.*res.multiplier, f);
% FFT_sf(isnan(FFT_sf)) = 0;
% 
% for w = 1:length(bw) 
%     center = win_center(w); 
%     win_location = (f >= center -(bw(w)/2)) & (f <= center + (bw(w)/2)); 
%     hann_window = hann(sum(win_location)); 
%     center_ind_f = dsearchn(f', center); 
%     segment = zeros(1, length(FFT_sf)); 
%     segment(1, center_ind_f-floor(length(hann_window)/2):center_ind_f+floor(length(hann_window)/2)-1) = hann_window'; 
%     segment = segment .* FFT_sf; 
%     
%     % First bin is f=0, then you have even number of bins
%     % Need to take first half and conjugate mirror to fill other.
%     Nhalf = (N-1)/2;
%     nonzeroHalf_sf =  segment(2:(2+Nhalf-1));
%     segment((Nhalf+2):end) = conj(nonzeroHalf_sf(end:-1:1));
%     
%     impulse_sf = ifftshift(ifft(segment));
%     
% %     nonzeroHalf_nf =  FFT_nf(2:(2+Nhalf-1));           %same for nf
% %     FFT_nf((Nhalf+2):end) = conj(nonzeroHalf_nf(end:-1:1));
% %     
% %     impulse_nf = ifftshift(ifft(FFT_nf));
%     
%     % Make time vector: t=0 will be at the center owing to ifftshift
%     t = (0:(N-1))/fs - timedur/2; % in milleseconds
%     
%     % Start a bit negative because D component has close to 0 delay and go to
%     % half of the reconstructed duration (50 ms) as planned (just for plotting)
%     t_min = -5e-3;
%     t_max = 15e-3;
%     inds_valid = t > t_min & t < t_max;
%     impulse_sf = impulse_sf(inds_valid);
%     t_sf = t(inds_valid);
%     
%     impulse_nf = impulse_nf(inds_valid);                %same for nf
%     
%     % Plot envelope of impulse response to see if there are two peaks, with the
%     % notch between peaks being somewhere in the 1-5 ms range.
%     figure(40);
%     impulse_sf_env = abs(hilbert(impulse_sf));
%     impulse_nf_env = abs(hilbert(impulse_nf));
%     plot(t_sf*1e3,  impulse_sf_env, 'linew', 2);
%     hold on;
%     plot(t_sf*1e3, impulse_nf_env, 'linew', 2)
%     xlim([t_min*1e3, 20]);
%     xlabel('Time (ms)', 'FontSize', 16);
%     ylabel('Impulse Response Envelope', 'FontSize', 16);
%     set(gca, 'FontSize', 16);
%     title('IFFT method',  'FontSize', 16);
%     
%     % Do windowing of signals
%     % % Start with  hard windows (box) and then smooth edges by 0.5 ms
%     % smoothing_kernel = blackman(ceil(0.5e-3*fs));
%     % smoothing_kernel  = smoothing_kernel / sum(smoothing_kernel);
%     
%     win_min = .5e-3; %seconds   Window to take as primary component, less reflections
%     win_max = 5e-3;
%     
%     t_sf_clean = t_sf.*(t_sf <= win_max & t_sf >= win_min);
%     impulse_sf_clean = impulse_sf.*(t_sf <= win_max & t_sf >= win_min);
%     impulse_nf_clean = impulse_nf.*(t_sf <= win_max & t_sf >= win_min);
%     
%     
%     
%     
% end




%%% HARI VERSION
complex_sf_raw = oae_complex .* res.multiplier; 
complex_nf_raw = noise_complex .* res.multiplier; 
fs = stim.Fs; % Exact value doesn't matter, but simplest to match orig

% Reconstruct time-domain response to 100 ms. Then the first 50 ms can be
% used, with the recognition that part of the right half is negative time.
% 50 ms should still be plenty for all OAE components to have come back
% and the impulse response to have decayed to noise floor.
timedur = 30e-3;

% Check if 1/(frequency resolution) is long enough
if 1/mean(abs(diff(testfreq))) < timedur/2
    warning('Too few DPOAE points, aliasing  will likely  occur');
end
N = roundodd(timedur * fs);
f = (0:(N-1))*fs/N; % FFT bin frequencies

% Window the frequency domain response to avoid sharp edges while filling
% in zeros for empty bins
rampsamps = 8;
ramp = hanning(2*rampsamps);
complex_sf_ramp  = complex_sf_raw;
complex_sf_ramp(1:rampsamps) = complex_sf_ramp(1:rampsamps)...
    .*ramp(1:rampsamps);
complex_sf_ramp((end-rampsamps+1):end) = ...
    complex_sf_ramp((end-rampsamps+1):end).* ramp((end-rampsamps+1):end);

complex_nf_ramp  = complex_nf_raw;              %same for nf
complex_nf_ramp(1:rampsamps) = complex_nf_ramp(1:rampsamps)...
    .*ramp(1:rampsamps);
complex_nf_ramp((end-rampsamps+1):end) = ...
    complex_nf_ramp((end-rampsamps+1):end).* ramp((end-rampsamps+1):end);

% Fill in FFT bins from dp data
FFT_sf =  interp1(testfreq, complex_sf_ramp, f);
FFT_sf(isnan(FFT_sf)) = 0;

FFT_nf =  interp1(testfreq, complex_nf_ramp, f); %same for nf
FFT_nf(isnan(FFT_nf)) = 0;

% First bin is f=0, then you have even number of bins
% Need to take first half and conjugate mirror to fill other.
Nhalf = (N-1)/2;
nonzeroHalf_sf =  FFT_sf(2:(2+Nhalf-1));
FFT_sf((Nhalf+2):end) = conj(nonzeroHalf_sf(end:-1:1));

impulse_sf = ifftshift(ifft(FFT_sf));

nonzeroHalf_nf =  FFT_nf(2:(2+Nhalf-1));           %same for nf
FFT_nf((Nhalf+2):end) = conj(nonzeroHalf_nf(end:-1:1));

impulse_nf = ifftshift(ifft(FFT_nf));

% Make time vector: t=0 will be at the center owing to ifftshift
t = (0:(N-1))/fs - timedur/2; % in milleseconds

% Start a bit negative because D component has close to 0 delay and go to
% half of the reconstructed duration (50 ms) as planned (just for plotting)
t_min = -5e-3;
t_max = 15e-3;
inds_valid = t > t_min & t < t_max;
impulse_sf = impulse_sf(inds_valid);
t_sf = t(inds_valid);

impulse_nf = impulse_nf(inds_valid);                %same for nf

% Plot envelope of impulse response to see if there are two peaks, with the
% notch between peaks being somewhere in the 1-5 ms range.
figure(40);
impulse_sf_env = abs(hilbert(impulse_sf));
impulse_nf_env = abs(hilbert(impulse_nf)); 
plot(t_sf*1e3,  impulse_sf_env, 'linew', 2);
hold on; 
plot(t_sf*1e3, impulse_nf_env, 'linew', 2)
xlim([t_min*1e3, 20]);
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Impulse Response Envelope', 'FontSize', 16);
set(gca, 'FontSize', 16);
title('IFFT method',  'FontSize', 16);

% Do windowing of signals
% % Start with  hard windows (box) and then smooth edges by 0.5 ms
% smoothing_kernel = blackman(ceil(0.5e-3*fs));
% smoothing_kernel  = smoothing_kernel / sum(smoothing_kernel);

win_min = .1e-3; %seconds   Window to take as primary component, less reflections
win_max = 21e-3; 

t_sf_clean = t_sf.*(t_sf <= win_max & t_sf >= win_min); 
impulse_sf_clean = impulse_sf.*(t_sf <= win_max & t_sf >= win_min); 
impulse_nf_clean = impulse_nf.*(t_sf <= win_max & t_sf >= win_min); 

% Error using conv2
% First and second arguments must be single or double.
% original code (which works elsewhere?) 
% win_D_only = conv(t_dp < D_only_dur, smoothing_kernel, 'same');
% win_R_only = conv(t_dp > D_only_dur, smoothing_kernel, 'same');

% trying the following to solve problem: 
% t_dp_D_only = t_sf.*(t_sf < D_only_dur); 
% t_dp_R_only = t_sf.*(t_sf > D_only_dur); 
% win_D_only = conv(t_dp_D_only, smoothing_kernel, 'same');
% win_R_only = conv(t_dp_R_only, smoothing_kernel, 'same');

% impulse_dp_D_only = impulse_sf .* win_D_only;
% impulse_dp_R_only = impulse_sf .* win_R_only;

complex_sf_allbins = fft(impulse_sf_clean);
f_complex_sf_allbins = (0:(numel(impulse_sf_clean)-1)) ...
    * fs/numel(impulse_sf_clean);
phasor_tmin_correction = ...
    exp(1j*2*pi*f_complex_sf_allbins*abs(t_min));
complex_sf_allbins_corrected = complex_sf_allbins ...
    .* phasor_tmin_correction;
complex_sf_IFFT = interp1(f_complex_sf_allbins,...
    complex_sf_allbins_corrected, testfreq);

complex_nf_allbins = fft(impulse_nf_clean);         %same for noise floor
f_complex_nf_allbins = (0:(numel(impulse_nf_clean)-1)) ...
    * fs/numel(impulse_nf_clean);
phasor_tmin_correction = ...
    exp(1j*2*pi*f_complex_nf_allbins*abs(t_min));
complex_nf_allbins_corrected = complex_nf_allbins ...
    .* phasor_tmin_correction;
complex_nf_IFFT = interp1(f_complex_nf_allbins,...
    complex_nf_allbins_corrected, testfreq);



%% Plot resulting figure (in SPL) // IFFT method
figure;
semilogx(testfreq/1000, db(abs(oae_complex).*res.multiplier), 'linew', 2);
hold on;
semilogx(testfreq/1000, db(abs(complex_sf_IFFT)), 'linew', 2);
hold on;
plot(testfreq/1000, db(abs(noise_complex).*res.multiplier), 'linew', 1.5);
plot(testfreq/1000, db(abs(complex_nf_IFFT)), '--', 'linew', 2)
title([subj, ' | SFOAE | ', condition], 'FontSize', 14 )
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
xlabel('Frequency (kHz)','FontWeight','bold')
ylabel('Amplitude (dB SPL)','FontWeight','bold')
legend('SFOAE', 'IFFTcleaned', 'NF', 'clean NF')







%% Qerb

% theta = unwrap(angle(oae_complex))/(2*pi); % cycles
% tau_pg = -diff(theta)./diff((testfreq')); %sec
% x = (testfreq(2:end) + testfreq(1:end-1))/2; 
% Nsf = x' .* (tau_pg);
% r = 1.25; % not exact
% Qerb = Nsf .* r; 

% Using IFFT cleaned version
theta = unwrap(angle(complex_sf_IFFT'))/(2*pi); % cycles
tau_pg = diff(theta)./diff((testfreq')); %sec
x = (testfreq(2:end) + testfreq(1:end-1))/2; 
Nsf = x' .* (tau_pg);
r = 1.25; % not exact
Qerb = Nsf .* r; 




% %% Plot resulting figure (in SPL)
% figure;
% plot(testfreq/1000, db(abs(oae_complex).*res.multiplier), 'linew', 1.75);
% hold on;
% plot(testfreq/1000, db(abs(noise_complex).*res.multiplier), '--', 'linew', 1.5);
% title([subj, ' | SFOAE | ', condition], 'FontSize', 14 )
% set(gca, 'XScale', 'log', 'FontSize', 14)
% xlim([.5, 16])
% xticks([.5, 1, 2, 4, 8, 16])
% xlabel('Frequency (kHz)','FontWeight','bold')
% ylabel('Amplitude (dB SPL)','FontWeight','bold')
% legend('SFOAE', 'NF')

%% Get EPL units
[SF] = calc_EPL(testfreq, complex_sf_IFFT', res.calib, 1);
% [SF] = calc_EPL(testfreq, oae_complex.*res.multiplier, res.calib, 1);
res.complex.sf_epl = SF.P_epl;
res.f_epl = SF.f;
res.dbEPL_sf = db(abs(SF.P_epl));

[NF] = calc_EPL(testfreq, complex_nf_IFFT', res.calib, 1);
% [NF] = calc_EPL(testfreq, noise_complex.*res.multiplier, res.calib, 1);
res.complex.nf_epl = NF.P_epl;
res.f = NF.f;
res.dbEPL_nf = db(abs(NF.P_epl));

%% Plot figure again (in EPL)
figure;
plot(testfreq/1000, res.dbEPL_sf, 'linew', 3, 'Color', 'k');
hold on;
plot(testfreq/1000, res.dbEPL_nf, 'k--', 'linew', 1.5);
title([subj ' | SFOAE | ' condition], 'FontSize', 14)
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 50])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('F2 Frequency (kHz)', 'FontWeight', 'bold')
legend('SFOAE', 'NF')
drawnow;

%% Summary Points of OAE amplitude
dpoae_full = res.dbEPL_sf;
dpnf_full = res.dbEPL_nf;
f2 = res.f/1000;

% SEt params
fmin = 0.5;
fmax = 16;
edges = 2 .^ linspace(log2(fmin), log2(fmax), 21);
bandEdges = edges(2:2:end-1);
centerFreqs = edges(3:2:end-2);

dpoae = zeros(length(centerFreqs),1);
dpnf = zeros(length(centerFreqs),1);
dpoae_w = zeros(length(centerFreqs),1);
dpnf_w = zeros(length(centerFreqs),1);
% resample / average to 9 center frequencies
for z = 1:length(centerFreqs)
    band = find( f2 >= bandEdges(z) & f2 < bandEdges(z+1));
    
   % NF from which SNR was calculated included median of 7 points
    % nearest the target frequency.
    for y = 1:length(band)
        dpnf_medians(y) = median(dpnf_full(band(y)-3:band(y)+3)); 
    end
        
    SNR = dpoae_full(band) - dpnf_medians;
    weight = (10.^(SNR./10)).^2;
    
    dpoae(z, 1) = mean(dpoae_full(band));
    dpnf(z,1) = mean(dpnf_full(band));
    
    dpoae_w(z,1) = sum(weight.*dpoae_full(band))/sum(weight);
    dpnf_w(z,1) = sum(weight.*dpnf_full(band))/sum(weight);
    
end


figure;
hold on;
semilogx(f2, dpoae_full, 'Color', [.8, .8, .8], 'linew', 2)
semilogx(f2, dpnf_full, '--', 'linew', 1.5, 'Color', [.8, .8, .8])
semilogx(centerFreqs, dpoae_w, 'o', 'linew', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 50])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('Frequency (kHz)', 'FontWeight', 'bold')

result.f = f2; 
result.oae_full = dpoae_full; 
result.nf_full = dpnf_full; 
result.centerFreqs = centerFreqs; 
result.oae_summary = dpoae_w; 
 
%% Save result function
res.windowdur = windowdur;
res.offsetwin = offsetwin;
res.npoints = npoints;
res.avgSFOAEresp = SFOAE;   % average mic response
res.t_freq = t_freq;
res.f = testfreq;           % frequency vectors
res.a = a;                  % coefficients
res.b = b;
res.tau = tau;
res.phasor = phasor;
res.subj = subj;
res.multiplier = stim.VoltageToPascal.* stim.PascalToLinearSPL;
res.complex.oae = oae_complex; %not IFFT cleaned
res.complex.nf = noise_complex;
res.durs = durs; 

res.theta = theta; 
res.tau_pg = tau_pg; 
res.f_pg = x; 
res.Nsf = Nsf; 
res.Qerb = Qerb; 
res.r = r; 

data.result = result; 
data.res = res;
%% Export:
if sedated_flag
    if ~exist([datapath, filesep, 'Processed', filesep, 'Sedated' filesep], 'dir')
        mkdir([datapath, filesep, 'Processed', filesep, 'Sedated' filesep])
    end
    cd([datapath, filesep, 'Processed', filesep, 'Sedated', filesep]);
else
    cd([datapath, filesep, 'Processed', filesep]);
end

fname = [subj,'_SFOAEswept_',condition, file(end-24:end-4) ];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'data')
cd(cwd);
