subj = input('Please subject ID:', 's');
% Make directory to save results
paraDir = './DEMO_Results/';
if(~exist(strcat(paraDir,'/',subj),'dir'))
    mkdir(strcat(paraDir,'/',subj));
end
useTDT = true;
rankList = [2];
%1 for left, 2 for right, 3 for both
ear = 2;
respDir = strcat(paraDir,'/',subj,'/');

%% Call app

%need to test this...
for i = 1:length(rankList)
    pause(1);
    rank = rankList(i)
    save('currsubj', 'subj', 'ear', 'respDir', 'useTDT','rank');
    f0dlAPP = DEMO_F0DL_HR_3AFC_2016b;
    while isvalid(f0dlAPP)
        pause(0.1);
    end
end

