% SFOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: October 27, 2023
% Purpose:
% Helpful info: Need to add Qerb calculation and consider eliminating long
% components (w/ IFFT method)

%%% Set these parameters %%%%%%%%%%%%%

windowdur = 0.038; % 40ms in paper
offsetwin = 0.00; % 20ms in paper
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
cwd = pwd;
cd(datapath)
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

stim = data.stim;

% SET CALIB FILE HERE
calib = data.FPL.FPLearData;
res.calib = calib; 

cd(cwd);

if ~isfield(stim, 'scale')
    stim.scale = 'log';
    stim.nearfreqs = [1.1,1.12, 1.14,1.16];
end

figure; plot(stim.SuppBuffs(1,1:400)); hold on; plot([128, 247], [0,0], 'or')
text(128, .1, '128'); text(247, .1, '247')
ask_delay = inputdlg('extra delay?');  % 247; %128
delay_oops = str2double(ask_delay{1}); 


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
if strcmp(stim.scale, 'log')    %linear sweep
    testfreq = 2 .^ linspace(log2(f1), log2(f2), npoints);
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else
    testfreq = linspace(f1, f2, npoints);
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end

%duration changes w/ frequency
durs = .038*(2.^(-0.3*(t_freq-stim.buffdur))-1)/ (-0.3*log(2)) + 0.038;

%% Artifact rejection

% Cancel out stimulus
SFOAEtrials = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs;
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

SFOAE = mean(resp_AR, 'omitNaN'); % mean SFOAE after artifact rejection

nfreqs = stim.nearfreqs; 
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



%% Qerb

theta = unwrap(angle(oae_complex))/(2*pi); % cycles
tau_pg = -diff(theta)./diff((testfreq')); %sec
x = (testfreq(2:end) + testfreq(1:end-1))/2; 
Nsf = x' .* (tau_pg);
r = 1.25; % not exact
Qerb = Nsf .* r; 




%% Plot resulting figure (in SPL)
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

%% Get EPL units
[SF] = calc_EPL(testfreq, oae_complex.*res.multiplier, res.calib, 1);
res.complex.sf_epl = SF.P_epl;
res.f_epl = SF.f;
res.dbEPL_sf = db(abs(SF.P_epl));

[NF] = calc_EPL(testfreq, noise_complex.*res.multiplier, res.calib, 1);
res.complex.nf_epl = NF.P_epl;
res.f = NF.f;
res.dbEPL_nf = db(abs(NF.P_epl));

%% Plot figure again (in EPL)
figure;
plot(testfreq/1000, res.dbEPL_sf, 'linew', 3, 'Color', '#4575b4');
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
    
    % Do some weighting by SNR
    
    % TO DO: NF from which SNR was calculated included median of 7 points
    % nearest the target frequency.
    SNR = dpoae_full(band) - dpnf_full(band);
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
semilogx(centerFreqs, dpoae_w, 'o', 'linew', 4, 'MarkerSize', 10, 'MarkerFaceColor', '#d73027', 'MarkerEdgeColor', '#d73027')
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
res.complex.oae = oae_complex;
res.complex.nf = noise_complex;
res.durs = durs; 

data.result = result; 
data.res = res;
%% Export:
cd(datapath);
fname = [subj,'_SFOAEswept_',condition, file(end-24:end-4) ];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'res')
cd(cwd);
