clear; clc; clf; close all
%load('D:\Animals\Earl\2021\07\20210706\E20210706.mat');

%How frequent were catch trials?
load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
numTrials = size(Data,2);
isCatch = 0;
for trial = 1:numTrials
    trialName = Data(trial).Overview.trialName;
    if strcmpi(trialName,'GridReachingCatch')
        isCatch = isCatch + 1;
    end
end
pctCatch = 100*isCatch/numTrials;

%Ever used "Center Hold Long" or "Target Hold Long" interval?
usedCHL = 0;
usedTHL = 0;
for trial = 1:numTrials
    centerHoldIntervalName = Data(trial).Parameters.StateTable(2).Interval.name{1,1};
    targetHoldIntervalName = Data(trial).Parameters.StateTable(5).Interval.name{1,1};
    if strcmpi(centerHoldIntervalName,'Center Hold Long')
        usedCHL = usedCHL + 1;
    end
    if strcmpi(targetHoldIntervalName,'Target Hold Long')
        usedTHL = usedTHL + 1;
    end
end

