info.measure = 'MST';
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
    cd Jane
end

info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;
subj = visit.subj.ID;

% Make directory to save results
paraDir = 'C:\Experiments\Sam\Jane\RESULTS\';
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');
save('currsubj', 'subj', 'respDir', 'info');

%% Call app
pause(1);
INFOMASK;