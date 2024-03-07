% SFOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: August 1, 2023
% Purpose:
% Helpful info:

%%% Set these parameters %%%%%%%%%%%%%

windowdur = 0.040; % 40ms in paper
offsetwin = 0.0; % 20ms in paper
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
cwd = pwd;
cd(datapath)
datafile = dir(fullfile(cd,['SFOAEswept_', subj, '*.mat']));
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

stim = data.stim;
resp = data.resp; 

%% Set variables needed from stim.
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
if abs(stim.speed) < 20 %linear sweep
    testfreq = 2 .^ linspace(log2(f1), log2(f2), npoints);
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else % log sweep
    testfreq = linspace(f1, f2, npoints);
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end

%duration changes w/ frequency
durs = .038*(2.^(-0.3*t_freq)-1)/ (-0.3*log(2)) + 0.038;

%% Artifact rejection

% Cancel out stimulus
SFOAEtrials = resp.ProbeBuffs + resp.SuppBuffs - resp.BothBuffs;
trials = size(SFOAEtrials,1);

% high pass filter the response (can also be done on ER10X hardware)
filtcutoff = 300;
b = fir1(1000, filtcutoff*2/stim.Fs, 'high');
SFOAEtrials= filtfilt(b, 1, SFOAEtrials')';

% Set empty matricies for next steps
coeffs = zeros(npoints, 6);
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
        if oae(j,k) > median_oae(1,k) + 4*std_oae(1,k)
            win = find( (t > (t_freq(k) - durs(k).*.1)) & ...
                (t < (t_freq(k) + durs(k).*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end


%% Calculate the Noise Floor (two ways)
numOfTrials = floor(trials/2)*2; % need even number of trials
% first method, subtraction
for y = 1:2:numOfTrials
    pos_SF(ceil(y/2),:) = resp_AR(y,:);
    neg_SF(ceil(y/2),:) = resp_AR(y+1,:);
end

numTrials2 = floor(numOfTrials/4).*2;
pos_noise = zeros(numTrials2/2, size(pos_SF,2));
neg_noise = zeros(numTrials2/2, size(pos_SF,2));
for x = 1:2:numTrials2
    pos_noise(ceil(x/2),:) = (pos_SF(x, :) - pos_SF(x+1, :)) / 2;
    neg_noise(ceil(x/2),:) = (neg_SF(x, :) - neg_SF(x+1, :)) / 2;
end

noise = [pos_noise; neg_noise];

SFOAE = mean(resp_AR, 'omitNaN'); % mean SFOAE after artifact rejection
NOISE = mean(noise, 'omitNaN'); % mean SFOAE after artifact rejection

%% LSF Analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
coeffs_n = zeros(npoints, 2);
tau = zeros(npoints, 1);
coeffs_noise = zeros(npoints,8);

durs = .038*(2.^(-0.3*t_freq)-1)/ (-0.3*log(2)) + 0.038;

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
    
    model_noise = [cos(1.1*phiProbe_inst(win)) .* taper;
        -sin(1.1*phiProbe_inst(win)) .* taper;
        cos(1.12*phiProbe_inst(win)) .* taper;
        -sin(1.12*phiProbe_inst(win)) .* taper;
        cos(1.14*phiProbe_inst(win)) .* taper;
        -sin(1.14*phiProbe_inst(win)) .* taper;
        cos(1.16*phiProbe_inst(win)) .* taper;
        -sin(1.16*phiProbe_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 2);
    coeff_n = zeros(maxoffset, 2);
    resid = zeros(maxoffset, 1);
    coeff_noise = zeros(maxoffset, 8);
    
    for offset = 0:maxoffset
        resp = SFOAE(win+offset) .* taper;
        resp_n = NOISE(win+offset) .* taper;
        
        coeff(offset + 1, :) = model' \ resp';
        coeff_n(offset + 1, :) = model' \ resp_n';
        coeff_noise(offset +1, :) = model_noise' \ resp';
        
        resid(offset +1) = sum( (resp - coeff(offset+1, :) * model).^2);
    end
    
    [~, ind] = min(resid);
    
    coeffs(k, :) = coeff(ind, :);
    coeffs_n(k, :) = coeff_n(ind, :);
    coeffs_noise(k,:) = coeff_noise(ind,:);
    
    tau(k) = (ind - 1) * (1/stim.Fs); % delay in sec
    
end

%% Amplitude and delay calculations
a = coeffs(:, 1);
b = coeffs(:, 2);
a_n = coeffs_n(:, 1); % subtraction nf
b_n = coeffs_n(:, 2);

phi = tau.*testfreq'; % cycles (from delay/offset)
phasor = exp(-1j * phi* 2 * pi);

% for noise
noise2 = zeros(npoints,4);
for i = 1:2:8
    noise2(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end

oae_complex = complex(a, b).*phasor;
noise_complex2 = complex(a_n, b_n);
noise_complex = mean(noise2,2);
res.multiplier = stim.VoltageToPascal.* stim.PascalToLinearSPL;



%% Windowing 
%% Separating D and R components by IFFT
% Code from IFFT method of separating D and R components from HB. 

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

win_min = 0; %seconds   Window to take as primary component, less reflections
win_max = 15e-3; 

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







%% Plot resulting figure
figure;
plot(testfreq/1000, db(abs(oae_complex).*res.multiplier), 'linew', 1.75);
hold on;
plot(testfreq/1000, db(abs(noise_complex).*res.multiplier), '--', 'linew', 1.5);
title([subj, ' | SFOAE | ', condition], 'FontSize', 14 )
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
xticks([.5, 1, 2, 4, 8, 16])
xlabel('Frequency (kHz)','FontWeight','bold')
ylabel('Amplitude (dB SPL)','FontWeight','bold')
legend('SFOAE', 'NF')

% %% Apply EPL

if subj(1) == 'S'
    calib1 = 0; 
    calib2 = 0; 
    % find calib file
    calibdir = [prefix '/THESIS/Pitch_Diagnostics_Data/FPLcalib/Human/' condition '/' subj];
    cd(calibdir)
    
    checkDIR=dir(sprintf('Calib_Ph1ER-10X_%s*.mat', subj));
    if isempty(checkDIR)
        fprintf('No such Calibs for subject %s\n',subj);
    elseif size(checkDIR, 1) > 1
        %fprintf('Date/Time of Test: %s\n',stim.date); 
        checkDIR =uigetfile(sprintf('Calib_Ph1ER-10X_%s*.mat', subj));
        load(checkDIR);
    else
        load(checkDIR.name);
    end
    
    res.calib.Ph1 = calib;
    calib1 = 1;
    clear checkDIR; 
    clear calib;

    
%     checkDIR=dir(sprintf('Calib_Ph2ER-10X_%s*%s*.mat', subj, stim.date(1:11)));
%     if isempty(checkDIR)
%         fprintf('No such Calibs for subject %s\n',subj);
%         return
%     elseif size(checkDIR, 1) > 1
%         fprintf('Date/Time of Test: %s\n',stim.date);
%         checkDIR =uigetfile(sprintf('Calib_Ph12ER-10X_%s*%s*.mat', subj, stim.date(1:11)));
%         load(checkDIR);
%     else
%         load(checkDIR.name);
%     end
%     
%     res.calib.Ph2 = calib;
%     calib2 = 1;
%     clear calib;
     cd(cwd)
    
    % if calibration is here
    if calib1 == 1 
        

        % Get EPL units
        [SF] = calc_EPL(testfreq, oae_complex.*res.multiplier, res.calib.Ph1);
        res.complex.sf_epl = SF.P_epl;
        res.f_epl = SF.f;
        res.dbEPL_sf = db(abs(SF.P_epl));
        
        [NF] = calc_EPL(testfreq, noise_complex.*res.multiplier, res.calib.Ph1);
        res.complex.nf_epl = NF.P_epl;
        res.f_epl = NF.f;
        res.dbEPL_nf = db(abs(NF.P_epl));
        
        % plot figure again
        figure;
        plot(testfreq/1000, res.dbEPL_sf, 'linew', 3, 'Color', '#d73027');
        hold on;
        plot(testfreq/1000, res.dbEPL_nf, 'k--', 'linew', 1.5);
        %plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
        %plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
        %title(sprintf('Subj: %s, Ear: %s', string(subj), string(ear)))
        title([subj ' | SFOAE | ' condition], 'FontSize', 14)
        set(gca, 'XScale', 'log', 'FontSize', 14)
        xlim([.5, 16])
        ylim([-50, 50])
        xticks([.5, 1, 2, 4, 8, 16])
        ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
        xlabel('F2 Frequency (kHz)', 'FontWeight', 'bold')
        legend('SFOAE', 'NF')
        drawnow;
    end
    
end

%
%         

%         
%         % plot figure again
%         figure;
%         plot(testfreq/1000, db(abs(SF.P_epl)));
%         hold on;
%         plot(testfreq/1000, db(abs(NF.P_epl)));
%         title(sprintf('SFOAE Subj: %s, Ear: %s', string(subj), string(ear)))
%         set(gca, 'XScale', 'log', 'FontSize', 14)
%         xlim([.5, 16])
%         xticks([.5, 1, 2, 4, 8, 16])
%         xlabel('Frequency (kHz)')
%         ylabel('Amplitude dB EPL')
%         legend('SFOAE', 'NF')
%     end
%     
% end

%% Save result function
res.windowdur = windowdur;
res.offsetwin = offsetwin;
res.npoints = npoints;
res.avgSFOAEresp = SFOAE;   % average mic response
res.avgNOISEresp = NOISE;
res.t_freq = t_freq;
res.f = testfreq;           % frequency vectors
res.a = a;                  % coefficients
res.b = b;
res.a_n = a_n;
res.b_n = b_n;
res.stim = stim;
res.tau = tau;
res.phasor = phasor;
res.subj = subj;
res.multiplier = stim.VoltageToPascal.* stim.PascalToLinearSPL;
res.complex.oae = oae_complex;
res.complex.nf = noise_complex;
res.complex.nf2 = noise_complex2;

%% Export:
cd(datapath);
fname = [subj,'_SFOAEswept_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'res')
cd(cwd);
