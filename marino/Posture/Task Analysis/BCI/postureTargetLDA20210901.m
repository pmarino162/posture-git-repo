clear; clc; clf; close all;

%% Load, Preprocess, and Label Data    
    [Data] = loadEarlData20210901;
    Data = Data([Data.trialStatus]==1);
    
%% Exlude trials that are too long
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
    [Data] = excludeLengths(Data,trialInclStates);

%% Run GPFA On States of Interest; Save projection to Spikes Field
    %Info for getStatesTraj    
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Center Target','Step 1'};
        trialInclStates(1).inclOccurrence = {'last','last'};
        trialInclStates(1).addTimeToBeginning = {-100,0};
        trialInclStates(1).addTimeToEnd = {0,100};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};      
        trialInclStates(2).inclStates = {'Center Hold','Delay','No Cheating','Target Acquire','Target Hold'};
        trialInclStates(2).inclOccurrence = {'last','last','last','last','last'};
        trialInclStates(2).addTimeToBeginning = {0,0,0,0,0};
        trialInclStates(2).addTimeToEnd = {0,0,0,0,0};
        

    %Info for GPFA
    numDims = 17;
    binWidth = 25;
    D = struct('trialId',[],'spikes',[]);
    numTrials = size(Data,2);
    D = repmat(D,1,numTrials);
    startTimes = zeros(1,numTrials);
    for trial = 1:numTrials
       D(trial).trialId = trial;
       [FR,FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSpikeBins','timeRelToTrialStart',true);
       D(trial).spikes = FR';
       startTimes(trial) = FRTimestamps(1);
    end
    % Run neuralTraj
    result = neuralTraj(001,D,'binWidth',binWidth,'xDim',numDims); 
    CSGPFAParams = result;
    clearvars D
    % Orthonormalize GPFA results and add to Data Struct
    trialIdList = [result.seqTrain.trialId];
    for trial = 1:numTrials
        trialInd = find(trialIdList == trial);
        [Xorth, Corth, TT]= orthogonalize(result.seqTrain(trialInd).xsm,result.estParams.C);
        GPFAProj = Xorth';
        GPFALength = result.seqTrain(trialInd).T;
        binStartTime = startTimes(trial) + binWidth/2;
        binEndTime = binStartTime+binWidth*(GPFALength-1);
        GPFABinTimes = binStartTime:binWidth:binEndTime;
        Data(trial).spikes.GPFA = GPFAProj;
        Data(trial).spikes.GPFABinTimes = GPFABinTimes;
    end
    CSGPFAParams.Corth = Corth;
    clearvars result
    clearvars allSpikeBins GPFABinTimes GPFALength GPFAProj trialIdList TT Xorth Corth trialInd
    rmdir mat_results s
    
%% Create trajStruct
    condFields = {{'task','conditionData','taskID'},{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'GPFA','marker'};
  
    %Delay
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'DelayedCenterOut20210828'};   
        trialInclStates(1).inclStates = {'Delay'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {0};    
    conditionData = [Data.conditionData]; 
    delayData = Data([conditionData.taskID]==2);
    kinData = [delayData.kinData];
    delayData = delayData([kinData.delayLength]>=500);
    delayTrajStruct = getTrajStruct(delayData,condFields,trajFields,trialInclStates);
    
    %Execution
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'BCI Center Out'};
        trialInclStates(1).inclStates = {'Step 1'};
        trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {100};
    trialInclStates(2).trialName = {'DelayedCenterOut20210828'};   
        trialInclStates(2).inclStates = {'No Cheating','Target Acquire'};
        trialInclStates(2).inclOccurrence = {'last','last'};
        trialInclStates(2).addTimeToBeginning = {-100,0};
        trialInclStates(2).addTimeToEnd = {0,100};    
        
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
            
%% Delay Target LDA
    task = 2;
    tempTrajStruct = delayTrajStruct([delayTrajStruct.task]==task);
    targetList = unique([tempTrajStruct.target]);
    numTargets = size(targetList,2);
    allTraj = []; allTargetLabels = [];
    for i = 1:size(tempTrajStruct,2)
       target = tempTrajStruct(i).target;
       numTraj = size(tempTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = tempTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
        end
    end
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    delayTargetLDA = targetLDA(:,1:numTargets-1);

%% Reach Target LDA
    task = 2;
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    targetList = unique([tempTrajStruct.target]);
    numTargets = size(targetList,2);
    allTraj = []; allTargetLabels = [];
    for i = 1:size(tempTrajStruct,2)
       target = tempTrajStruct(i).target;
       numTraj = size(tempTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = tempTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
        end
    end
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    reachTargetLDA = targetLDA(:,1:numTargets-1);

%% BC Target & Posture LDA
    task = 1;
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = unique([tempTrajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    allTraj = []; allPostureLabels = []; allTargetLabels = [];
    for i = 1:size(tempTrajStruct,2)
       target = tempTrajStruct(i).target;
       posture = tempTrajStruct(i).posture;
       numTraj = size(tempTrajStruct(i).allGPFA,2);
       for j = 1:numTraj
            %Trajectory
            traj = tempTrajStruct(i).allGPFA(j).traj;
            numSteps = size(traj,1);
            allTraj = vertcat(allTraj,traj);
            %Labels
            trajTargetLabel = ones(numSteps,1).*target;
            allTargetLabels = vertcat(allTargetLabels,trajTargetLabel);
            trajPostureLabel = ones(numSteps,1).*posture;
            allPostureLabels = vertcat(allPostureLabels,trajPostureLabel);
        end
    end
    
    %Elbow Posture and Target
    elbowPostureInds = allPostureLabels==4 | allPostureLabels == 5;
    elbowPostureLabels = allPostureLabels(elbowPostureInds);
    elbowTargetLabels = allTargetLabels(elbowPostureInds);
    elbowTraj = allTraj(elbowPostureInds,:);
    
    elbowPostureLDA = fisherLDA(elbowTraj, elbowPostureLabels);
    [elbowPostureLDA,~] = qr(elbowPostureLDA);
    elbowPostureLDA = elbowPostureLDA(:,1);
    
    %Shoulder Posture and Target
    shoulderPostureInds = allPostureLabels==1 | allPostureLabels == 3;
    shoulderPostureLabels = allPostureLabels(shoulderPostureInds);
    shoulderTargetLabels = allTargetLabels(shoulderPostureInds);
    shoulderTraj = allTraj(shoulderPostureInds,:);
    
    shoulderPostureLDA = fisherLDA(shoulderTraj, shoulderPostureLabels);
    [shoulderPostureLDA,~] = qr(shoulderPostureLDA);
    shoulderPostureLDA = shoulderPostureLDA(:,1);
    
    %All Posture and Target
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    BCPostureLDA = postureLDA(:,1:numPostures-1);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    BCTargetLDA = targetLDA(:,1:numTargets-1);
    
    [postTargOrth,~] = qr([BCPostureLDA(:,1:3),BCTargetLDA(:,1:2)]);
    BCPostTargOrth = postTargOrth(:,1:5);
    

%% Add LDA Projections to Traj Struct and Model Struct
    %Delay
    for i = 1:size(delayTrajStruct,2)
       delayTrajStruct(i).avgDelayTargetLDA = delayTrajStruct(i).avgGPFA.traj*delayTargetLDA;
       delayTrajStruct(i).avgReachTargetLDA = delayTrajStruct(i).avgGPFA.traj*reachTargetLDA;
       delayTrajStruct(i).avgBCPostureLDA = delayTrajStruct(i).avgGPFA.traj*BCPostureLDA;
       delayTrajStruct(i).avgBCTargetLDA = delayTrajStruct(i).avgGPFA.traj*BCTargetLDA;
       delayTrajStruct(i).avgBCPostTargOrth = delayTrajStruct(i).avgGPFA.traj*BCPostTargOrth;
       delayTrajStruct(i).shoulderPostureLDA = delayTrajStruct(i).avgGPFA.traj*shoulderPostureLDA;
       delayTrajStruct(i).elbowPostureLDA = delayTrajStruct(i).avgGPFA.traj*elbowPostureLDA;
    end
    %Reach
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgDelayTargetLDA = trajStruct(i).avgGPFA.traj*delayTargetLDA;
       trajStruct(i).avgReachTargetLDA = trajStruct(i).avgGPFA.traj*reachTargetLDA;
       trajStruct(i).avgBCPostureLDA = trajStruct(i).avgGPFA.traj*BCPostureLDA;
       trajStruct(i).avgBCTargetLDA = trajStruct(i).avgGPFA.traj*BCTargetLDA;
       trajStruct(i).avgBCPostTargOrth = trajStruct(i).avgGPFA.traj*BCPostTargOrth;
       trajStruct(i).shoulderPostureLDA = trajStruct(i).avgGPFA.traj*shoulderPostureLDA;
       trajStruct(i).elbowPostureLDA = trajStruct(i).avgGPFA.traj*elbowPostureLDA;
    end

%% Load Colormaps
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    scmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customBlue.mat');
    ecmap = customBlue;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;

%% Plot Delay Neural and Marker Trajectories
    targetList = unique([delayTrajStruct.target]);

    %LDA Spaces
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).avgDelayTargetLDA;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).avgReachTargetLDA;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-2;
               case {5,6}
                 	traj = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.traj;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'Delay'; 'Targ LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'Delay'; 'Targ LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'Reach'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'Reach'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Delay')
    
    
    
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).shoulderPostureLDA;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {2}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).elbowPostureLDA;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-1;
               case {3,4}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).avgBCPostTargOrth;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.traj;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Shoulder LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Elbow LDA 1'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Delay')

    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).avgBCPostTargOrth;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = delayTrajStruct([delayTrajStruct.target]==target).avgBCPostTargOrth;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.traj;
                    time = delayTrajStruct([delayTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Post LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Post LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Delay')
    
%% Plot Reach Neural and Marker Trajectories
    reachTrajStruct = trajStruct([trajStruct.task]==2);
    targetList = unique([reachTrajStruct.target]);

    %LDA Spaces
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).avgDelayTargetLDA;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).avgReachTargetLDA;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-2;
               case {5,6}
                 	traj = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.traj;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'Delay'; 'Targ LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'Delay'; 'Targ LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'Reach'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'Reach'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Reach')
    
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).shoulderPostureLDA;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {2}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).elbowPostureLDA;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-1;
               case {3,4}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).avgBCPostTargOrth;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.traj;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Shoulder LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Elbow LDA 1'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Reach')
    
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).avgBCPostTargOrth;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = reachTrajStruct([reachTrajStruct.target]==target).avgBCPostTargOrth;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.traj;
                    time = reachTrajStruct([reachTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Post LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Post LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-85 85])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Reach')
    
%% Plot BC Neural and Marker Trajectories
    BCTrajStruct = trajStruct([trajStruct.task]==1 & [trajStruct.posture]==3);
    targetList = unique([BCTrajStruct.target]);
    hcCenterLoc = [-25,-420];
    %LDA Spaces
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).avgDelayTargetLDA;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).avgReachTargetLDA;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-2;
               case {5,6}
                 	traj = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.traj-hcCenterLoc;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'Delay'; 'Targ LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'Delay'; 'Targ LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'Reach'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'Reach'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-215 -45])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Brain Control')
    
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).shoulderPostureLDA;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {2}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).elbowPostureLDA;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim-1;
               case {3,4}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).avgBCPostTargOrth;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.traj-hcCenterLoc;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Shoulder LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Elbow LDA 1'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                ylim([-215 -45])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Brain Control')
    
    
    f = figure; f.Position = [10 10 400 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.07,0.2);
    hold on
    for target = targetList
        for dim = 1:6
           axes(ha(dim));
           switch dim
               case {1,2}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).avgBCPostTargOrth;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim;
               case {3,4}
                    traj = BCTrajStruct([BCTrajStruct.target]==target).avgBCPostTargOrth;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgGPFA.timestamps;
                    dimInd = dim+1;
               case {5,6}
                 	traj = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.traj-hcCenterLoc;
                    time = BCTrajStruct([BCTrajStruct.target]==target).avgMarker.timestamps;
                    dimInd = dim-4;
           end
           plot(time,traj(:,dimInd),'LineWidth',2,'Color',tcmap(target,:));
           hold on
           switch dim
            case 1
                ylabel({'BC'; 'Post LDA 1'})
                ylim([-1 1])
            case 2
                ylabel({'BC'; 'Post LDA 2'})
                ylim([-1 1])
            case 3
                ylabel({'BC'; 'Targ LDA 1'})
                ylim([-1 1])
            case 4
                ylabel({'BC'; 'Targ LDA 2'})
                ylim([-1 1])
            case 5
                ylabel({'Hand Pos'; 'x (mm)'})
                 ylim([-215 -45])
            case 6
                ylabel({'Hand Pos'; 'y (mm)'})
                ylim([-85 85])
                xlabel('time (ms)')
          end
        end
    end
    sgtitle('Brain Control')
    
%% Angle Bt Axes 
dot(shoulderPostureLDA,elbowPostureLDA)


%% 3D Views of BCI Trajectories
	task  = 1;     
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = unique([tempTrajStruct.target]);
    
    figure
    xDim = 1; yDim = 2; zDim = 3;
    for posture = postureList
        switch posture
            case 1
                cmap = scmap(1,:);
            case 3
                cmap = scmap(5,:);
            case 4
                cmap = ecmap(1,:);
            case 5
                cmap = ecmap(5,:);
        end
        for target = targetList
            postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
            plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(1,:),'LineWidth',2); 
            hold on
            plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
        end
    end
    xlabel(['Posture LDA 1']); ylabel(['Posture LDA 2']); zlabel(['Target LDA 1'])
    grid on
    axis equal

%%
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Decomposition Analysis\Earl20200317_Decomposition Analysis\';
    
    for task = 1:2      
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        postureList = unique([tempTrajStruct.posture]);
        targetList = unique([tempTrajStruct.target]);
        
        %Top 3 GPFA - Actual
        figure
        xDim = 1; yDim =2; zDim = 3;
        for posture = postureList
            for target = targetList
                traj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.traj;
                plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
            end
        end
        xlabel(['GPFA ',num2str(xDim)]); ylabel(['GPFA ',num2str(yDim)]); zlabel(['GPFA ',num2str(zDim)])
        grid on
        titleStr = ['Top 3 GPFA']; title(titleStr);
        %saveas(gcf,[dirStr,titleStr,'.fig'])

        %Posture/Target Orth Projection - Actual
        figure
        xDim = 1; yDim = 3; zDim = 4;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal
        %saveas(gcf,[dirStr,titleStr,'.fig'])     

        figure
        xDim = 2; yDim = 3; zDim = 4;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA 2']); ylabel(['Target LDA 1']); zlabel(['Target LDA 2'])
        grid on
        axis equal
        
        figure
        xDim = 1; yDim = 2; zDim = 3;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Posture LDA 2']); zlabel(['Target LDA 1'])
        grid on
        axis equal
        
        figure
        xDim = 1; yDim = 2; zDim = 3;
        for posture = postureList
            for target = targetList
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                plot3(postTargOrthProj(1:end-1,xDim),postTargOrthProj(1:end-1,yDim),postTargOrthProj(1:end-1,zDim),'Color',cmap(posture,:),'LineWidth',2); 
                hold on
                plot3(postTargOrthProj(1,xDim),postTargOrthProj(1,yDim),postTargOrthProj(1,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(5,xDim),postTargOrthProj(5,yDim),postTargOrthProj(5,zDim),'.','Color',tcmap(target,:),'MarkerSize',20); 
                plot3(postTargOrthProj(10,xDim),postTargOrthProj(10,yDim),postTargOrthProj(10,zDim),'.','Color',tcmap(target,:),'MarkerSize',20);  
            end
        end
        xlabel(['Posture LDA 1']); ylabel(['Posture LDA 2']); zlabel(['Target LDA 1'])
        grid on
        axis equal
%         %PxT Projection - Actual
%         figure
%         xDim = 1; yDim = 2; zDim = 3;
%         for posture = postureList
%             for target = targetList
%                 pxtProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgPxTLDA;
%                 plot3(pxtProj(:,xDim),pxtProj(:,yDim),pxtProj(:,zDim),'Color',cmap(posture,:),'LineWidth',2); 
%                 hold on
%             end
%         end
%         xlabel(['PxT LDA ',num2str(xDim)]); ylabel(['PxT LDA ',num2str(yDim)]); zlabel(['PxT LDA ',num2str(zDim)])
%         grid on
%         axis equal
%         titleStr = ['PxT']; title(titleStr);
    
    end

    
    %% Posture, Target, and Null vs Time

    for task = 1:2     
        f = figure; f.Position = [10 10 900 300];
        [ha, pos] = tight_subplot(2,5,0.15,0.15,0.05);
        tempTrajStruct = trajStruct([trajStruct.task]==task);
        postureList = unique([tempTrajStruct.posture]);
        targetList = [1,5];
    
        for posture = postureList
            targetInd = 0;
            for target = targetList
                targetInd = targetInd + 1;
                postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
                nullProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgNull;
                time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
                for sp = 1:4
                   axes(ha(sp));
                   hold on
                   if targetInd == 1
                        plot(time,postTargOrthProj(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
                   elseif targetInd == 2
                       plot(time,postTargOrthProj(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
                   end
                end
                for sp = 5:10
                    axes(ha(sp));
                   hold on
                   if targetInd == 1
                        plot(time,nullProj(:,sp-4),'Color',cmap(posture,:),'LineWidth',2);   
                   elseif targetInd == 2
                       plot(time,nullProj(:,sp-4),':','Color',cmap(posture,:),'LineWidth',2);
                   end
                end
            end
        end
        sgtitle(['Task ',num2str(task)])
        for sp = 1:10
            axes(ha(sp)); ax = gca; ax.TickDir = 'out';
            xticks([-100 0 250]); xticklabels({'-100',0,'250'})
            ylimits = ylim; ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
%             ylimits = [-2.5, 2.5];
            ax.YLim = ylimits;
            yticks([ylimits(1) ylimits(2)]); yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
            xlabel('time after go cue (ms)')
            switch sp
                case 1
                    ylabel(['Posture LDA 1 (a.u.)'])
                case 2
                    ylabel(['Posture LDA 2 (a.u.)'])
                case 3
                    ylabel(['Target LDA 1 (a.u.)'])
                case 4
                    ylabel(['Target LDA 2 (a.u.)'])
               
                otherwise
                    ylabel(['Null ',num2str(sp-3),' (a.u.)'])
            end
            set(gca,'fontname','arial')
        end
% 
%     
%         % GPFA vs time
%         f = figure; f.Position = [10 10 1000 400];
%         [ha, pos] = tight_subplot(2,6,0.1,0.1,0.05);
%         for posture = postureList
%             for target = targetList
%                 GPFA = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.traj;
%                 time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
%                 for sp = 1:12
%                    axes(ha(sp));
%                    hold on
%                    if target == targetList(1)
%                         plot(time,GPFA(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
%                    elseif target == targetList(2)
%                        plot(time,GPFA(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
%                    end
%                 end
%             end
%         end

    end
    
%% Compare time courses to marker time courses 
%     condFields = {{'task','conditionData','taskID'},{'posture','conditionData','postureID'},{'target','targetData','targetID'}};
%     trajFields = {'GPFA','marker'};
%     
%     conditionData = [Data.conditionData];
%     taskID = [conditionData.taskID];
%     task2Data = Data(taskID==2);
%     
%  
%     trialInclStates(2).trialName = {'Posture Device Active Movements'};
%         trialInclStates(2).inclStates = {'No Cheating','Target Acquire','Target Hold'};
%         trialInclStates(2).inclOccurrence = {'last','last','last'};
%         trialInclStates(2).addTimeToBeginning = {0,0,0};
%         trialInclStates(2).addTimeToEnd = {0,0,0};    
%         
%     markerTrajStruct = getTrajStruct(task2Data,condFields,trajFields,trialInclStates);
% 




    task = 2;     
    f = figure; f.Position = [10 10 500 700];
    [ha, pos] = tight_subplot(6,1,0.05,0.05,0.15);
    tempTrajStruct = trajStruct([trajStruct.task]==task);
    postureList = unique([tempTrajStruct.posture]);
    targetList = [1,5];
    
    for posture = 2
        targetInd = 0;
        for target = targetList
            targetInd = targetInd + 1;
            postTargOrthProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgBCPostTargOrth;
            nullProj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgNull;
            time = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgGPFA.timestamps;
            markerTraj = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgMarker.traj;
            markerTime = tempTrajStruct(find([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target)).avgMarker.timestamps;
            for sp = 1:4
               axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(time,postTargOrthProj(:,sp),'Color',cmap(posture,:),'LineWidth',2); 
               elseif targetInd == 2
                   plot(time,postTargOrthProj(:,sp),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
            for sp = 5:6
               axes(ha(sp));
               hold on
               if targetInd == 1
                    plot(markerTime,markerTraj(:,sp-4),'Color',cmap(posture,:),'LineWidth',2);   
               elseif targetInd == 2
                   plot(markerTime,markerTraj(:,sp-4),':','Color',cmap(posture,:),'LineWidth',2);
               end
            end
        end
    end
    sgtitle(['Task ',num2str(task)])
        
        
    for sp = 1:6
        axes(ha(sp)); ax = gca; ax.TickDir = 'out';
        xticks([-100 0 250]); xticklabels({'-100',0,'250'})
        ylimits = ylim; ylimits = [round(ylimits(1),1), round(ylimits(2),1)];
%             ylimits = [-2.5, 2.5];
        ax.YLim = ylimits;
        yticks([ylimits(1) ylimits(2)]); yticklabels({num2str(ylimits(1)), num2str(ylimits(2))})
        xlabel('time after go cue (ms)')
        switch sp
            case 1
                ylabel(['Posture LDA 1 (a.u.)'])
            case 2
                ylabel(['Posture LDA 2 (a.u.)'])
            case 3
                ylabel(['Target LDA 1 (a.u.)'])
            case 4
                ylabel(['Target LDA 2 (a.u.)'])

            otherwise
                ylabel(['Null ',num2str(sp-3),' (a.u.)'])
        end
        set(gca,'fontname','arial')
    end
% 
