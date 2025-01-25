clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220401 - time-varying posture analysis 1';
        
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;

%% Load Data, get trajStruct
    [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForTrajDist(Data,'reaching');
    
    trajFields = {'smoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;
    trialInclStates(1).trialName = {'GridReaching'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    %Hold
        trialInclStates(1).inclStates = {{'state','Center Hold','first',50},{'state','Delay','first',0}};
        holdTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false,'getTrialAverages',true);
    %Plan 
        kinData = [Data.kinData];
        delayLength = [kinData.delayLength];
        rmTrials = delayLength < 500;
        trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Delay','first',250}};
        planTrajStruct = getTrajStruct20211210(Data(~rmTrials),condFields,trajFields,trialInclStates,binWidth,'matchConditions',false,'getTrialAverages',true);
    %Reach
        trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
        reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false,'getTrialAverages',true);
    %%Whole Trial
        %trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Success with Reward','first',0}};
        %wholeTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
%% Get posture and target lists
    tempTrajStruct = holdTrajStruct;
        postureList = unique([tempTrajStruct.posture]); 
        numPostures = size(postureList,2);
        targetList = unique([tempTrajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(tempTrajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(tempTrajStruct,2);
    
%% Do posture LDA on on all phases
    for phase = {'hold','plan','reach','all'}
        switch phase{1,1}
            case 'hold'
                phaseTrajStruct = holdTrajStruct;
            case 'plan'
                phaseTrajStruct = planTrajStruct;
            case 'reach'
                phaseTrajStruct = reachTrajStruct;
            case 'all'
                phaseTrajStruct = [holdTrajStruct,planTrajStruct];
        end
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        structInd = 1;   
        for posture = postureList
            posture
           obsStruct(structInd).label = posture;
           allObs = [];
           tempTrajStruct = phaseTrajStruct([phaseTrajStruct.posture]==posture);
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   traj = tempTrajStruct(i).allSmoothFR(j).trialAvg;
                   %traj = tempTrajStruct(i).allSmoothFR(j).traj;
                   allObs = vertcat(allObs,traj);
               end
           end
           obsStruct(structInd).allObs = allObs;
           obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
           structInd = structInd + 1;
        end
        
        switch phase{1,1}
            case 'hold'
                [holdPostureLDA] = doLDA(obsStruct);
            case 'plan'
                [planPostureLDA] = doLDA(obsStruct);
            case 'reach'
                [reachPostureLDA] = doLDA(obsStruct);
            case 'all'
                [allPostureLDA] = doLDA(obsStruct);
        end
    end
    
%% Create result struct
    pts = struct('posture',[],'pts',[]);
    ellipse = struct('posture',[],'mu',[],'ellipse',[]);
    resultStruct = struct('trainPhase','','testPhase','','pts',pts,'ellipse',ellipse);
    structInd = 1;
    for trainPhase = {'hold','plan','reach','all'}
        switch trainPhase{1,1}
            case 'hold'
               proj = holdPostureLDA;
            case 'plan'
               proj = planPostureLDA;
            case 'reach'
               proj = reachPostureLDA;
            case 'all'
                proj = allPostureLDA;
        end
        for testPhase = {'hold','plan','reach','all'}
            switch testPhase{1,1}
               case 'hold'
                   trajStruct = holdTrajStruct;
               case 'plan'
                   trajStruct = planTrajStruct;
               case 'reach'
                   trajStruct = reachTrajStruct;
               case 'all'
                   trajStruct = [holdTrajStruct,planTrajStruct];  
            end 
            postureInd = 1;
            pts = struct('posture',[],'pts',[]);
            ellipse = struct('posture',[],'mu',[],'a',[]);
            for posture = postureList
                traj = [trajStruct([trajStruct.posture]==posture).allSmoothFR];
                traj = vertcat(traj.trialAvg);
                ldaProj = traj*proj;
                ldaProj = ldaProj(:,1:2);
                pts(postureInd).posture = posture;
                pts(postureInd).pts = ldaProj;
                ellipse(postureInd).posture = posture;
                ellipse(postureInd).mu = mean(ldaProj);
                ellipse(postureInd).ellipse = error_ellipse(ldaProj);
                postureInd = postureInd + 1;
            end
            resultStruct(structInd).trainPhase = trainPhase{1,1};
            resultStruct(structInd).testPhase = testPhase{1,1};
            resultStruct(structInd).pts = pts;
            resultStruct(structInd).ellipse = ellipse;
            structInd = structInd + 1;
        end
    end
    
% %% Visualize All Data
%     for trainPhase = {'hold','plan','reach'}
%         for testPhase = {'hold','plan','reach'}
%             figure
%             for posture = postureList
%                 
%                switch testPhase{1,1}
%                    case 'hold'
%                        trajStruct = holdTrajStruct;
%                    case 'plan'
%                        trajStruct = planTrajStruct;
%                    case 'reach'
%                        trajStruct = reachTrajStruct;
%                end
%                traj = [trajStruct([trajStruct.posture]==posture).allSmoothFR];
%                traj = vertcat(traj.trialAvg);
%                switch trainPhase{1,1}
%                    case 'hold'
%                        proj = holdPostureLDA;
%                    case 'plan'
%                        proj = planPostureLDA;
%                    case 'reach'
%                        proj = reachPostureLDA;
%                end
%                ldaProj = traj*proj;
%                plot3(ldaProj(:,1),ldaProj(:,2),ldaProj(:,3),'.','MarkerSize',20,'Color',pcmap(posture,:));
%                hold on
%             end
%             title(['Train ',trainPhase{1,1},'Test ',testPhase{1,1}])
%             xlabel('LDA 1')
%             ylabel('LDA 2')
%             zlabel('LDA 3')
%         end
%     end
%     
%% 2d plot
    f = figure;
    f.Position = [0 0 800 800]
    plotInd = 1;
    for trainPhase = {'hold','plan','reach','all'}
        for testPhase = {'hold','plan','reach','all'}
            subplot(4,4,plotInd)
            ellipse = resultStruct(strcmpi({resultStruct.trainPhase},trainPhase{1,1}) & strcmpi({resultStruct.testPhase},testPhase{1,1})).ellipse;
            for posture = postureList
               mu = ellipse([ellipse.posture]==posture).mu;
               ell = ellipse([ellipse.posture]==posture).ellipse;
               plot(mu(:,1),mu(:,2),'.','MarkerSize',20,'Color',pcmap(posture,:));
               hold on
               plot(mu(:,1)+ell(:,1),mu(:,2)+ell(:,2),'-','Color',pcmap(posture,:));
            end
            title(['Train ',trainPhase{1,1},' Test ',testPhase{1,1}])
            %xlabel('LDA 1')
            %ylabel('LDA 2')
            plotInd = plotInd + 1;
        end
    end


%% 2d reach trajectory plot 
    figure
    for posture = [2]
    for target = [1:8]
        traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj;
        traj = traj*allPostureLDA;
        traj = traj(4:20,:);
        plot3(traj(:,1),traj(:,2),traj(:,3),'LineWidth',3,'Color',tcmap(target,:))
        hold on
        plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',20,'Color',tcmap(target,:))
        plot3(traj(end,1),traj(end,2),traj(end,3),'>','Markersize',7,'Color',tcmap(target,:))
    end
    end
    xlabel('LDA 1')
    ylabel('LDA 2')
    zlabel('LDA 3')
    axis equal
%% Local function for performing LDA
     function [LDAproj] = doLDA(obsStruct)
        %Preallocate
        minNumObs = min([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 
        %Fill
        k = 1;
        for i = 1:size(obsStruct,2)
            totalClassObs = size(obsStruct(i).allObs,1);
            obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
            labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
            k = k+minNumObs;
        end
        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end  