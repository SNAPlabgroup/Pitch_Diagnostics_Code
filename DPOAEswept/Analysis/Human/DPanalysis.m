% DPOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: December 11, 2023
% Purpose:
% Helpful info:

%%%%%%%%% Set these parameters %%%%%%%%%%%%%%%%%%

windowdur = .25; 
offsetwin = 0.0; % not finding additional delay
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
cwd = pwd;
cd(datapath)
datafile = dir(fullfile(cd,['DPOAEswept_', subj, '_*.mat']));
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
if exist('stim')
    stim.scale = 'log';
    stim.nearfreqs = [0.9,.88, .86,.84];
else
    stim = data.stim; 
    resp = data.resp; 
end

%% Get appropriate calibration file

% find calib file
calibdir = [prefix, 'THESIS', filesep,...
    'Pitch_Diagnostics_Data', filesep, 'FPLcalib', filesep, 'Human', ...
    filesep, condition, filesep, subj];
cd(calibdir)

checkDIR=dir(sprintf('Calib_Ph1ER-10X_%s*%s*.mat', subj, data.info.date(1:11)));
if isempty(checkDIR)
    fprintf('No such Calibs for subject %s\n',subj);
elseif size(checkDIR, 1) > 1
    fprintf('Date/Time of Test: %s\n',data.info.date);
    checkDIR =uigetfile(sprintf('Calib_Ph1ER-10X_%s*%s*.mat', subj, data.info.date(1:11)));
    load(checkDIR);
else
    load(checkDIR.name);
end

calib.Ph1 = calib;

cd(cwd)

%% Analysis set up

% Set variables from the stim
phi1_inst = 2 * pi * stim.phi1_inst;
phi2_inst = 2 * pi * stim.phi2_inst;
phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
rdp = 2 / stim.ratio - 1;    % f_dp = f2 * rdp

trials = size(data.resp.AllBuffs,1);

t = stim.t;
if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else
    f_start = stim.fmin;
    f_end = stim.fmax;
end

% set freq we're testing and the timepoints when they happen.
if strcmp(stim.scale, 'log')        % in octave scaling
    freq_f2 = 2 .^ linspace(log2(f_start), log2(f_end), npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = log2(freq_f2/f_start)/stim.speed + stim.buffdur;
else                            % otherwise linear scaling
    freq_f2 = linspace(f_start, f_end, npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = (freq_f2-f_start)/stim.speed + stim.buffdur;
end

nfreqs = stim.nearfreqs;

%% Artifact Rejection

% high pass filter the response (can also be done on ER10X hardware)
%filtcutoff = 300;
%b = fir1(1000, filtcutoff*2/stim.Fs, 'high');
%DPOAEtrials= filtfilt(b, 1, stim.resp')';
DPOAEtrials = resp.AllBuffs;

% Set empty matricies for next steps
coeffs = zeros(npoints, 2);
a_temp = zeros(trials, npoints);
b_temp = zeros(trials, npoints);

% Least Squares fit of DP Only for AR
for x = 1:trials
    DPOAE = DPOAEtrials(x, :);
    fprintf(1, 'Checking trial %d / %d for artifact\n', x, (trials));
    
    for k = 1:npoints
        win = find( (t > (t_freq(k) - windowdur/2)) & ...
            (t < (t_freq(k) + windowdur/2)));
        taper = hanning(numel(win))';
        
        model_dp = [cos(phi_dp_inst(win)) .* taper;
            -sin(phi_dp_inst(win)) .* taper];
        
        resp = DPOAE(win) .* taper;
        
        coeffs(k, 1:2) = model_dp' \ resp';
    end
    a_temp(x,:) = coeffs(:, 1);
    b_temp(x,:) = coeffs(:, 2);
end

oae = abs(complex(a_temp, b_temp));
median_oae = median(oae);
std_oae = std(oae);

resp_AR = DPOAEtrials;
for j = 1:trials
    for k = 1:npoints
        if oae(j,k) > median_oae(1,k) + 3*std_oae(1,k)
            win = find( (t > (t_freq(k) - windowdur.*.1)) & ...
                (t < (t_freq(k) + windowdur.*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end

DPOAE = mean(resp_AR, 1, 'omitNaN');

%% LSF analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
tau_dp = zeros(npoints, 1); % delay if offset > 0
coeffs_noise = zeros(npoints,8);

% Least Squares fit of Chirp model (stimuli, DP, noise)
for k = 1:npoints
    
    fprintf(1, 'Running window %d / %d\n', k, npoints);
    
    % if using durs: windowdur = durs(k);
    win = find( (t > (t_freq(k) - windowdur/2)) & ...
        (t < (t_freq(k) + windowdur/2)));
    taper = hanning(numel(win))';
    
    % set the response
    resp = DPOAE(win) .* taper;
    
    % DP Coeffs with variable delay calculation
    model_dp = [cos(phi_dp_inst(win)) .* taper;
        -sin(phi_dp_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 6);
    coeff_n = zeros(maxoffset, 6);
    resid = zeros(maxoffset, 3);
    for offset = 0:maxoffset
        resp = DPOAE(win+offset) .* taper;
        coeff(offset + 1, 1:2) = model_dp' \ resp';
        resid(offset + 1, 1) = sum( (resp  - coeff(offset + 1, 1:2) * model_dp).^2);
    end
    [~, ind] = min(resid(:,1));
    coeffs(k, 1:2) = coeff(ind, 1:2);
    
    % Calculate delay
    tau_dp(k) = (ind(1) - 1) * 1/stim.Fs; % delay in sec
    
    % F1 Coeffs
    model_f1 = [cos(phi1_inst(win)) .* taper;
        -sin(phi1_inst(win)) .* taper];
    coeffs(k, 3:4) = model_f1' \ resp';
    
    % F2 Coeffs
    model_f2 = [cos(phi2_inst(win)) .* taper;
        -sin(phi2_inst(win)) .* taper];
    coeffs(k, 5:6) = model_f2' \ resp';
    
    % Noise Coeffs
    model_noise = ...
        [cos(nfreqs(1)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(1)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(2)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(2)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(3)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(3)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(4)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(4)*phi_dp_inst(win)) .* taper];
    coeffs_noise(k,:) = model_noise' \ resp';
    
end


%% Amplitude and Delay calculations
a_dp = coeffs(:, 1);
b_dp = coeffs(:, 2);
a_f1 = coeffs(:, 3);
b_f1 = coeffs(:, 4);
a_f2 = coeffs(:, 5);
b_f2 = coeffs(:, 6);

% complex DPOAE
oae_complex = complex(a_dp, b_dp);

% complex average noise
noise = zeros(npoints,4);
for i = 1:2:8
    noise(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end
noise_complex = mean(noise,2);

% delay
phi_dp = tau_dp.*freq_dp'; % cycles (from delay/offset)
phasor_dp = exp(-1j * phi_dp * 2 * pi);

VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;
res.VtoSPL = VtoSPL;

% %% Plot Results Figure
% figure;
% plot(freq_f2/1000, db(abs(oae_complex).*VtoSPL), 'linew', 2, 'Color', [0 0.4470 0.7410]);
% hold on;
% plot(freq_f2/1000, db(abs(noise_complex).*VtoSPL), '--', 'linew', 2, 'Color', [0.6350 0.0780 0.1840]);
% plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*VtoSPL), 'linew', 2, 'Color', [0.4940 0.1840 0.5560]);
% plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*VtoSPL), 'linew', 2, 'Color', [0.9290 0.6940 0.1250]);
% title('DPOAE', 'FontSize', 14)
% set(gca, 'XScale', 'log', 'FontSize', 14)
% xlim([.5, 16])
% ylim([-50, 90])
% xticks([.5, 1, 2, 4, 8, 16])
% ylabel('Amplitude (dB SPL)', 'FontWeight', 'bold')
% xlabel('F2 Frequency (kHz)', 'FontWeight', 'bold')
% legend('OAE', 'NF', 'F2', 'F1')
% drawnow;

%% Get EPL units
[DP] = calc_EPL(freq_dp, oae_complex.*VtoSPL, calib);
res.complex.dp_epl = DP.P_epl;
res.f_epl = DP.f;
res.dbEPL_dp = db(abs(DP.P_epl));

[NF] = calc_EPL(freq_dp, noise_complex.*VtoSPL, calib);
res.complex.nf_epl = NF.P_epl;
res.f_epl = NF.f;
res.dbEPL_nf = db(abs(NF.P_epl));

%% Summary Points of OAE amplitude
res.f.f2 = freq_f2;         % frequency vectors
res.f.f1 = freq_f1;
res.f.dp = freq_dp;

dpoae_full = res.dbEPL_dp;
dpnf_full = res.dbEPL_nf;
f2 = res.f.f2/1000;

%% Calculate Summary Points
% Set params
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
    weight = (10.^(SNR./10)).^2; % weighted averaged based on SNR
    
    dpoae(z, 1) = mean(dpoae_full(band));
    dpnf(z,1) = mean(dpnf_full(band));
    
    dpoae_w(z,1) = sum(weight.*dpoae_full(band))/sum(weight);
    dpnf_w(z,1) = sum(weight.*dpnf_full(band))/sum(weight);
    
end

%% Plot final figure 
figure;
hold on;
semilogx(f2, dpoae_full, 'Color', [.8, .8, .8], 'linew', 2)
semilogx(f2, dpnf_full, '--', 'linew', 1.5, 'Color', [.8, .8, .8])
semilogx(centerFreqs, dpoae_w, 'o', 'linew', 2, 'MarkerSize', 7, 'MarkerEdgeColor', 'b')
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 50])
xticks([.5, 1, 2, 4, 8, 16])
yticks([-50:10:50])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('F2 Frequency (kHz)', 'FontWeight', 'bold')
title(sprintf('%s | DPOAE | %s', subj, condition), 'FontSize', 16); 

%% Save final resulting variables
result.f2 = f2;
result.oae_full = dpoae_full;
result.nf_full = dpnf_full;
result.centerFreqs = centerFreqs;
result.oae_summary = dpoae_w;

res.windowdur = windowdur;
res.offsetwin = offsetwin;
res.npoints = npoints;
res.avgDPOAEresp = DPOAE;   % average mic response
res.t_freq = t_freq;
res.f.f2 = freq_f2;         % frequency vectors
res.f.f1 = freq_f1;
res.f.dp = freq_dp;
res.a.dp = a_dp;            % coefficients
res.b.dp = b_dp;
res.a.f1 = a_f1;
res.b.f1 = b_f1;
res.a.f2 = a_f2;
res.b.f2 = b_f2;
res.tau.dp = tau_dp;
res.complex.oae = oae_complex;
res.complex.nf = noise_complex;

data.result = result;
data.res = res;
data.stim = stim;
data.calib = calib; 
%% Export:
cd(datapath);
fname = [subj,'_DPOAEswept_',condition, file(end-24:end-4)];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'data')
cd(cwd);
