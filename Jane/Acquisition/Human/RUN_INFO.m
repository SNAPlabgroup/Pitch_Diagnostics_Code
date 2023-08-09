subj = input('Please subject ID:', 's');
% Make directory to save results
paraDir = 'C:\Experiments\Sam\SIN_INFO_Sam\RESULTS\';
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');
save('currsubj', 'subj', 'respDir');

%% Call app
pause(1);
INFOMASK;