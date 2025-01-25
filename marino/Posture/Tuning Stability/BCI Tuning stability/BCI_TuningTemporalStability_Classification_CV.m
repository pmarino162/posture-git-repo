clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, subselect
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    task = 'BCI';
    [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForTrajDist(Data,task);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
    trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numDims = size(trajStruct(1).avgSmoothFR.traj,2);
    timeList = 1:4;
    numTimes = size(timeList,2);
    
%% Look at histogram of trial lengths
    trialLengths = [];
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allSmoothFR,2)
            trialLengths = [trialLengths,size(trajStruct(i).allSmoothFR(j).timestamps,2)];
        end
    end
    histogram(trialLengths);
         
%% Preallocate obsStruct
    obsStruct = struct('target',[],'posture',[],'time',[],'postureLabelVec',[],'targetLabelVec',[],'timeLabelVec',[],'numObs',[],'allObs',[]);
    structInd = 1;
    for target = targetList
        for posture = postureList
            for time = timeList
                obsStruct(structInd).target = target;
                obsStruct(structInd).posture = posture;
                obsStruct(structInd).time = time;
                structInd = structInd + 1;
            end
        end
    end
    
%% Fill obsStruct
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allSmoothFR,2)
            timestamps = trajStruct(i).allSmoothFR(j).timestamps;
            trialLength = timestamps(end);
            if trialLength >= 200 && trialLength <=650
                target = trajStruct(i).target;
                posture = trajStruct(i).posture;
                traj = trajStruct(i).allSmoothFR(j).traj;
                for time = timeList
                    startTime = (time-1)*100;
                    endTime = time*100 - 25;
                    startInd = find(timestamps==startTime);
                    endInd = find(timestamps==endTime);
                    if ~isempty(startInd) && ~isempty(endInd)
                        trialAvg = mean(traj(startInd:endInd,:),1);
                        structInd = find([obsStruct.target]==target & [obsStruct.posture]==posture & [obsStruct.time]==time);
                        obsStruct(structInd).allObs = vertcat(obsStruct(structInd).allObs,trialAvg);
                    end
                end
            end
        end
    end
    for structInd = 1:size(obsStruct,2)
        numObs = size(obsStruct(structInd).allObs,1);
        obsStruct(structInd).numObs = numObs;
        obsStruct(structInd).postureLabelVec = ones(numObs,1)*obsStruct(structInd).posture;
        obsStruct(structInd).targetLabelVec = ones(numObs,1)*obsStruct(structInd).target;
        obsStruct(structInd).timeLabelVec = ones(numObs,1)*obsStruct(structInd).time;
    end

%% Get time-specific train and test sets for posture classification
    timeObsStruct = struct('time',[],'obs',[],'labels',[],'partition',[]);
    numFolds = 5;
    structInd = 1;
    for time = timeList
        minNumObs = min([obsStruct([obsStruct.time]==time).numObs]);
        obs = NaN(minNumObs*numPostures*numTargets,numDims);
        labels = NaN(minNumObs*numPostures*numTargets,1);
        i = 1;
        for target = targetList
            for posture = postureList
                numObs = obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture & [obsStruct.time]==time).numObs;
                obs(i:i+minNumObs-1,:) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture & [obsStruct.time]==time).allObs(randsample(numObs,minNumObs),:));
                labels(i:i+minNumObs-1,1) = vertcat(obsStruct([obsStruct.target]==target & [obsStruct.posture]==posture & [obsStruct.time]==time).postureLabelVec(1:minNumObs));
                i = i + minNumObs;
            end
        end
        c = cvpartition(minNumObs*numPostures*numTargets,'KFold',numFolds);
        timeObsStruct(structInd).time = time;
        timeObsStruct(structInd).obs = obs;
        timeObsStruct(structInd).labels = labels;
        timeObsStruct(structInd).partition = c;
        structInd = structInd + 1;
    end
   
%% Train posture classifier on individual time; classify across times
    acrossTimeResult = zeros(numTimes,numTimes,numFolds);
    for trainTime = timeList
       trainObs = timeObsStruct([timeObsStruct.time]==trainTime).obs;
       trainLabels = timeObsStruct([timeObsStruct.time]==trainTime).labels;
       trainPartition = timeObsStruct([timeObsStruct.time]==trainTime).partition;
       for testTime = timeList
            testObs = timeObsStruct([timeObsStruct.time]==testTime).obs;
            testLabels = timeObsStruct([timeObsStruct.time]==testTime).labels;
            testPartition = timeObsStruct([timeObsStruct.time]==testTime).partition;
            for fold = 1:numFolds
               [NBModel,numClasses] = NBTrainModel(trainObs(trainPartition.training(fold),:),trainLabels(trainPartition.training(fold),:)); 
               [predictedLabels,pctCorrect] = NBClassify(testObs(testPartition.test(fold),:),testLabels(testPartition.test(fold),:),numClasses,NBModel);
               acrossTimeResult(trainTime,testTime,fold) = pctCorrect;
            end   
       end
    end    
    acrossTimeResult = mean(acrossTimeResult,3);


    
%% Plot across-target posture classification performance    
    f = figure
    f.Position = [100 100 650 500];
    clims = [20,max(max(acrossTimeResult))];
    colormap hot
    imagesc(acrossTimeResult,clims)
    xticks([1:numTimes])
    xticklabels({'0-100','100-200','200-300','300-400'})
    xtickangle(45)
    xlabel('Testing Time (ms)')
    yticks([1:numTimes])
    yticklabels({'0-100','100-200','200-300','300-400'})
    ylabel('Training Time (ms)')
    c = colorbar;
    c.Label.String = 'Classification Accuracy (%)';
    ax = gca;
    ax.FontSize = 16;
