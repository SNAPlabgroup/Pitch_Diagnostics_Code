subj = input('Please subject ID:', 's');
% Make directory to save results
paraDir = './Results/';
if(~exist(strcat(paraDir,'/',subj),'dir'))
    mkdir(strcat(paraDir,'/',subj));
end
useTDT = true;
respDir = strcat(paraDir,'/',subj,'/');


% Additional parameters for display
infotext = 'In which interval did you hear a beep?';
prepad = 250e-3;
intervaldur = 250e-3;
isi = 400e-3;
postpad = 0;


save('currsubj', 'subj', 'respDir', 'useTDT',...
    'infotext', 'prepad', 'intervaldur', 'isi', 'postpad');

%% Call app
pause(1);
ToneThreshold3AFC_2016b;