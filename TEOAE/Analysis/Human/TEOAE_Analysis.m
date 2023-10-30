%% TEOAE Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: August 2, 2023
% Purpose:
% Helpful info:



%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,['CEOAE*',subj,'*.mat'])).name};

numOfFiles = length(datafile);
if numOfFiles > 1
    fprintf('More than one file...Select the right one?\n');
    datafile = uigetfile(['CEOAE*',subj,'*.mat']);
elseif numOfFiles < 1
    fprintf('No files for this subject...Quitting.\n')
    cd(cwd);
    return
else
    datafile = datafile{1};
end

load(datafile);

cd(cwd);

% %% Apply EPL


% find calib file
calibdir = ['/Volumes/SNH/THESIS/Pitch_Diagnostics_Data/FPLcalib/Human/' condition '/' subj];
cd(calibdir)

checkDIR=dir(sprintf('Calib_Ph1ER-10X_%s*.mat', subj));
if isempty(checkDIR)
    fprintf('No such Calibs for subject %s\n',subj);
elseif size(checkDIR, 1) > 1
    fprintf('Date/Time of Test: %s\n',click.date);
    checkDIR =uigetfile(sprintf('Calib_Ph1ER-10X_%s*.mat', subj));
    load(checkDIR);
else
    load(checkDIR.name);
end

res.calib.Ph1 = calib;
calib1 = 1;
clear checkDIR;
clear calib;


% Get EPL units
[TE] = calc_EPL(click.freq, click.Resp, res.calib.Ph1);
res.complex.teoae_epl = TE.P_epl;
res.f_epl = TE.f;
res.dbEPL_teoae = db(abs(TE.P_epl));

[NF] = calc_EPL(click.freq, click.NoiseFloor, res.calib.Ph1);
res.complex.nf_epl = NF.P_epl;
res.f_epl = NF.f;
res.dbEPL_nf = db(abs(NF.P_epl));

% plot figure again
figure;
plot(click.freq/1000, res.dbEPL_teoae, 'linew', 3, 'Color', '#d73027');
hold on;
plot(click.freq/1000, res.dbEPL_nf, 'k--', 'linew', 1.5);
%plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
%plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
%title(sprintf('Subj: %s, Ear: %s', string(subj), string(ear)))
title([subj ' | TEOAE | ' condition], 'FontSize', 14)
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-60, 0])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB EPL)', 'FontWeight', 'bold')
xlabel('Frequency (kHz)', 'FontWeight', 'bold')
legend('TEOAE', 'NF')
drawnow;


%% Analysis loop
% figure(1);
% hold on;
% plot(click.freq*1e-3, db(abs(click.Resp)), 'linew', 2);
% ylabel('Response (dB SPL)', 'FontSize', 14, 'FontWeight', 'bold');
% uplim = max(db(abs(click.Resp)));
% hold on;
% semilogx(click.freq*1e-3, db(abs(click.NoiseFloor)), 'linew', 2);
% xlabel('Frequency (kHz)', 'FontSize', 14, 'FontWeight', 'bold');
% legend('TEOAE', 'NoiseFloor');
% xlim([0.4, 16]);
% ticks = [0.5, 1, 2, 4, 8, 16];
% set(gca, 'XTick', ticks, 'FontSize', 14, 'xscale', 'log');
% ylim([-60, uplim + 5]);
% title([subj ' | TEOAE | ' condition], 'FontSize', 14)
% drawnow;

res.freq = click.freq;
res.db_resp = db(abs(click.Resp));
res.db_nf = db(abs(click.NoiseFloor));
%% Save Variables and figure
cd(datapath);
fname = [subj,'_TEOAE_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'res')
cd(cwd);

