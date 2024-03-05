% Analysis
%Take FFRs, process them, and calculate the cummulative sum of Harmonic
%Magnitudes. Should be generalized enough to run with any trial.
%Originally from: Andrew SivaprakaRAM, 03/2021
%Last modified: RAMantha Hauser, 03/2023


%% Parameters:

% ENTER
F0 = 223; %Fundamental freq of interest in Hz
harmonics = 8;

%% Import data
cwd = pwd;

cd(datapath)
datafile = dir(fullfile(cd,['*/p*.mat']));
cd(datafile(1).folder)
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

fname_out = [file(1:end-4),'_matlab.mat'];
cd(cwd);

%% Parameters
Fs0 = data.Stimuli.RPsamprate_Hz; 
Fs = 8e3;
iterations = 100;
window = [0.1,0.9];
gain = 20e3; %Gain (doesn't really matter since looking at SNR & PLV, used for magnitude accuracy)

K_MRS = 200; %number of distributions to average spectra and PLVs over
NF_iters = 100;

%% Analysis 

RAM_tot_full = data.AD_Data.AD_All_V{1, 1};
dbRAM = data.Stimuli.calib_dBSPLout - data.Stimuli.atten_dB;
trials = length(RAM_tot_full)/2;

%% Calculate the DFT for Responses
[RAM_f,RAM_DFT,RAM_PLV, floor_RAM, pos, neg] = getDFT(RAM_tot_full,trials,window,Fs,Fs0,gain,K_MRS,NF_iters);

for i = 1:length(pos)
    pos_waveforms(i,:) = pos{1,i};
    neg_waveforms(i,:) = neg{1,i};
end
res.pos_waveforms_all = pos_waveforms;
res.neg_waveforms_all = neg_waveforms;

res.pos_mean = mean(pos_waveforms, 1);
res.neg_mean = mean(pos_waveforms, 1);
res.env_waveform = res.pos_mean + res.neg_mean / 2;
res.tfs_waveform = res.pos_mean - res.neg_mean / 2;

res.f = RAM_f;
res.DFT = RAM_DFT;
res.PLV = RAM_PLV;
res.floor = floor_RAM;
% Plot of just Noise Floor
%         figure;
%         plot(RAM_f,floor_RAM);
%         title(strcat('Noise Floor - RAM:  ', f0str, ' Hz - Chin: ', chins(c)));
%         xlabel('Frequency');
%         ylabel('Attenuation');

%% Converting back to linear

RAM_DFT_uv = 10.^(RAM_DFT/20);


%% Plotting & Summation

        MAG = figure;
        subplot(2,1,1)
        hold on;
        plot(RAM_f,RAM_DFT)
        title(strcat('DFT with Noise Floor removed - RAM: 223 Hz - Chin: ', subj))
        %ylabel('SNR (dB)/Magnitude (dB, arbitrary)')
        ylabel('SNR (Linear Scale)')
        xlabel('Frequency')
        xlim([0,4e3])
        ylim([0,max(RAM_DFT*1.25)])


        %PLV Figure:
        PLV = figure;
        subplot(2,1,1)
        hold on;
        plot(RAM_f,RAM_PLV)
        title(strcat('PLV of Multiple Conditions - RAM: 223 Hz - Chin: ', subj))
        %ylabel('SNR (dB)/Magnitude (dB, arbitrary)')
        ylabel('PLV')
        xlabel('Frequency')
        xlim([0,2e3])
        ylim([0,1])

%Get peaks and sum them, look at crossings
% peaks with DFT method
[RAM_SUM,RAM_PKS,RAM_LOCS] = getSum(RAM_f,RAM_DFT_uv,F0,harmonics);

% peaks with PLV method
[SAMP_SUM,SAMP_PKS,SAMP_LOCS] = getSum(RAM_f,RAM_PLV,F0,harmonics);

RAM_MAG_SUM = RAM_SUM(end);

RAM_PLV_SUM = SAMP_SUM(end);

RAM_PLV_Rat = sum(SAMP_PKS(3:end))/sum(SAMP_PKS(1:2));


% save results to struct
res.DFT_uv = RAM_DFT_uv;
res.SUM_DFT = RAM_SUM;
res.PKS_DFT = RAM_PKS;
res.LOCS_DFT = RAM_LOCS;
res.SUM_PLV = SAMP_SUM;
res.PKS_PLV = SAMP_PKS;
res.LOCS_PLV = SAMP_LOCS;
res.MAG_SUM = RAM_MAG_SUM;
res.PLV_SUM = RAM_PLV_SUM;
res.PLV_Ratio = RAM_PLV_Rat;
res.MAG_MEAN = RAM_MAG_MEAN;
res.MAG_std = RAM_MAG_std;
res.PLV_MEAN = RAM_PLV_MEAN;
res.PLV_std = RAM_PLV_std;
res.PLV_MEAN_ratio = RAM_PLV_MEAN_rat;
res.PLV_std_ratio = RAM_PLV_std_rat;
res.chin = chin;
res.F0 = F0;
res.harms = harmonics;
res.iterations = iterations;
res.k = k;
res.window = window;
res.stim_dB = dbRAM;
res.Fs = Fs;
res.Fs0 = Fs0;


%% Saving Data

% Save individual chin data
savename = sprintf('Results_Q%s_%s.mat',chin, f0str);
cd(datadir)
save(savename, 'res');
fprintf('Results saved!\n');



