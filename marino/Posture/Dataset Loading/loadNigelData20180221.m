function [Data] = loadNigelData20180221

Data = [];
        %Posture 3
            prevTrialNum = 0;
            for blockInd = 54306:54325
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 3;
                    tempData(i).postureData.posture = 'A90';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Posture 2
            prevTrialNum = 0;
            for blockInd = 54331:54350
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 2;
                    tempData(i).postureData.posture = 'I45';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Posture 1
            prevTrialNum = 0;
            for blockInd = [54356:54369,54371:54375]
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 1;
                    tempData(i).postureData.posture = 'N00';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Posture 2
            prevTrialNum = 0;
            for blockInd = 54381:54400
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 2;
                    tempData(i).postureData.posture = 'I45';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Posture 3
            prevTrialNum = 0;
            for blockInd = 54406:54425
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 3;
                    tempData(i).postureData.posture = 'A90';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Posture 1
            prevTrialNum = 0;
            for blockInd = 54431:54525
                %Load/Preprocess Data Block
                loadStr = ['D:\Animals\Nigel\2018\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
                rawData = load(loadStr);
                [tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
                %Update trialNum; Add Posture
                for i=1:size(tempData,2)
                    tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                    tempData(i).postureData.postureID = 1;
                    tempData(i).postureData.posture = 'N00';
                end
                %Add to Data Struct; Update trialNum
                Data = [Data,tempData];
                prevTrialNum = Data(end).trialNum;
            end
        %Keep only Successful trials
            Data = Data([Data.trialStatus]==1);



%% Keep only trials with length w/in 2 stdDevs of mean length
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
    trialInclStates(1).inclStates = {'Cursor Release','Center Exit','Target Hold'};
    trialInclStates(1).inclOccurrence = {'last','last','last'};

    [Data] = excludeLengths(Data,trialInclStates);    
    
    
%% Label Conditions
    numTrials = size(Data,2);
    for trial = 1:numTrials
       Data(trial).Condition = ['P',num2str(Data(trial).postureData.posture),'T',num2str(Data(trial).targetData.targetID)]; 
    end

    
%% Get minimum number of successful trials for any condition
    conditionList = unique({Data.Condition});
    minNumber = size(Data,2);
    maxNumber = 0;
    for condition = conditionList
       tempData = Data(strcmpi({Data.Condition},condition));
       trialStatus = [tempData.trialStatus];
       numSuc = sum(trialStatus);
       if numSuc < minNumber
           minNumber = numSuc;
       end
       if numSuc > maxNumber
           maxNumber = numSuc;
       end
    end
    
    
%% Keep only minimum number
    minSucData = [];
    for condition = conditionList
        tempData = Data(strcmpi({Data.Condition},condition));
        tempData = tempData([tempData.trialStatus]==1);
        minSucData = [minSucData,tempData(1:minNumber)];
    end
    Data = minSucData;
    clearvars minSucData






end


