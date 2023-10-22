info.measure = 'MRT_tonecomplex';
info.version = 'v01';

if exist('C:\Experiments\Sam\current_visit.mat','file')
    load('C:\Experiments\Sam\current_visit.mat', 'visit')
    ask = questdlg(sprintf('Is this subject %s?', visit.subj.ID), ...
        'Check Subject', 'Yes', 'No', 'No');
else
    ask = 'No';
end

if strcmp(ask, 'No')
    cd ..
    startVisit
    cd(info.measure)
end

info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;
subj = visit.subj.ID;

% Make directory to save results
paraDir = 'C:\Experiments\Sam\MRT_tonecomplex\RESULTS\';
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end

startflag = questdlg('Start from the beginning?', 'Start', 'Yes', 'No', 'Yes');
if strcmp(startflag,'Yes')
    startBlock = 1;
else
    block = inputdlg('Enter the block number at which to start (1-13):', 'Start Block');
    startBlock = double(string(block{1})); 
end
respDir = strcat(paraDir,'\',subj,'\');
save('currsubj', 'subj', 'respDir', 'startBlock', 'info');

%% Call app
pause(1);
WiN6AFC;