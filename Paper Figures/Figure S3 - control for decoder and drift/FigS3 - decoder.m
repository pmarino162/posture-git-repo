clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S3 - control for decoder';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10;
    
%% Load data; get trajStruct
    dataset = 'E20190729';
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);  
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 

%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(numPostures);    

%% Fit posture PC on data with matched posture and decoder
    %Get posturePCs on blocks 1 and 2
    matchedTrajStruct = trajStruct([trajStruct.posture]==[trajStruct.decoder]); 
    %Project all data down to top PCs
    [allTraj,Mu] = collectAllAvgTraj(matchedTrajStruct);
    [matchedTrajStruct,PCs] = projectToTopPCs(allTraj,matchedTrajStruct,numPCsToKeep);
    % Fit P and T spaces, components for each group  
    [pSig,~] = getPandTsig(matchedTrajStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
    pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numPCsToKeep]);
    [pDims,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 

%% Add inidvidual trial projections to trajStruct
     for i = 1:size(trajStruct,2)
        numObs = size(trajStruct(i).allZSmoothFR,2);
        pcProj = zeros(numObs,2);
        for j = 1:numObs
           obs = mean(trajStruct(i).allZSmoothFR(j).traj); 
           pcProj(j,:) = (obs-Mu)*PCs(:,1:10)*pDims(:,1:2);
        end
        trajStruct(i).pcProj = pcProj;
     end

%% Plot and save
    fs = 14;
    %pcmap = pcmap([1,5],:);
    figure; hold on;
    plotInd = 1;
    for posture = 1:2
        for decoder = 1:2
            tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.decoder]==decoder);
            plotData = vertcat(tempTrajStruct.pcProj);
            jitter = 0.5*(rand(1,length(plotData))-.5);
            plot(plotInd*ones(1,length(plotData))+jitter,plotData(:,1),'.','MarkerSize',10,'Color',pcmap(posture,:));
            plotInd = plotInd + 1;
        end
    end
    %ylabel('Projection of neural activity \n onto Posture Dim 1 (a.u.)')
    xlim([0.5 4.5])
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ax= gca;
    ax.TickDir = 'out';
    %Save
    if saveFig
        saveas(gcf,fullfile(saveDir,[dataset,'_postureAxisVTrial.svg']));
    end


