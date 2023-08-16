subj = input('Please subject ID:', 's');
% Make directory to save results
paraDir = './Results/';
if(~exist(strcat(paraDir,'/',subj),'dir'))
    mkdir(strcat(paraDir,'/',subj));
end
useTDT = true;
randomize = true;
% rankList = [5,6,7,10,13,16];
rankList = 2:2:12;
rankList = [2:2:6,10,12];

%1 for left, 2 for right, 3 for both
ear = 2;
respDir = strcat(paraDir,'/',subj,'/');

%% Call app

if randomize
    rankList = rankList(randperm(length(rankList),length(rankList))); 
end

%need to test this...
for i = 1:length(rankList)
    pause(1);
    rank = rankList(i)
    save('currsubj', 'subj', 'ear', 'respDir', 'useTDT','rank');
    f0dlAPP = F0DL_HR_3AFC_Adaptive_2016b;
    while isvalid(f0dlAPP)
        pause(0.1);
    end
end


%% TODO:
% update newer version of Matlab script with transposed currentSignal
% update level in 2016 and new version
% needed to change the gui scaleable thing, after converting script to 2016
% compatible version
% check backwards compatibility...