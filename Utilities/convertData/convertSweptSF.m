function [filename, data] = convertSweptSF(x, subj, file, condition,gender)

info.measure = x.General.program_name(1:end-2);
info.version = x.General.version;
info.room = 'LYLE3035';
info.univ = 'Purdue';
info.researcher = file(1:2);
info.dir = file; 
info.date = x.General.date;
info.subj.ID = subj;
info.subj.status = condition;
info.gender = gender;
info.ear = x.sweptSFOAEData.stim.ear;
info.FPL_file = x.invfilterdata.coefFileNum; 
stim = x.sweptSFOAEData.stim;
resp.trialsCollected = x.sweptSFOAEData.stim.NUMtrials_Completed;
resp.ProbeBuffs = x.sweptSFOAEData.stim.ProbeBuffs;
resp.SuppBuffs = x.sweptSFOAEData.stim.SuppBuffs;
resp.BothBuffs = x.sweptSFOAEData.stim.BothBuffs;

data.info = info;
data.stim = stim;
data.resp = resp;

p_cal = dir(sprintf('p%04.f_calib_FPL_raw.mat', info.FPL_file));
y = load(p_cal.name, 'x'); 
data.FPL = y.x; 

time = x.General.time;
time(3:3:6) = '-';
filename = [info.measure '_' subj '_' info.date '_' time];

end