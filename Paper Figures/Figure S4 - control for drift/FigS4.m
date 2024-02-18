clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S4 - control for drift';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    numPCsToKeep = 10;    
    
%% Run loop for each dataset    
    for datasetList = {'N20171215','E20190830'}
        %Load data; get trajStruct
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);  
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
        % Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);    
        
        %Fit posture PCs
        %Project all data down to top PCs
        [allTraj,Mu] = collectAllAvgTraj(trajStruct);
        [trajStruct,PCs] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        % Fit P and T spaces, components for each group  
        [pSig,~] = getPandTsig(trajStruct,'avgPCA',minNumTimestamps,postureList,numPostures,targetList,numTargets);
        pSigReshape = reshape(squeeze(pSig),[numPostures*minNumTimestamps,numPCsToKeep]);
        [pDims,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 

        %Create trialStruct with individual trial observations
        trialStruct = struct('trial',[],'posture',[],'PC',[],'PCLDA',[]);
        structInd = 1;
        for i = 1:size(trajStruct,2)
            posture = trajStruct(i).posture;
            for j = 1:size(trajStruct(i).allZSmoothFR,2)
                trial = trajStruct(i).allZSmoothFR(j).trialNum;
                traj = trajStruct(i).allPCA(j).traj;
                traj = mean(traj,1);
                trialStruct(structInd).trial = trial;
                trialStruct(structInd).posture = posture;
                trialStruct(structInd).PC = traj*pDims(:,1);
                structInd = structInd + 1;
            end
        end
        
        %Plot and save
        fs = 14;
        %Individual trial observations in posture PC space
        figure; hold on;
        for i = 1:numel(trialStruct)
           trial = trialStruct(i).trial;
           posture = trialStruct(i).posture;
           PC = trialStruct(i).PC;
           plot(trial,PC,'.','MarkerSize',10,'Color',pcmap(posture,:));
        end
        %Add posture means for E20190830 dataset
        if strcmpi(dataset,'E20190830')
            for posture = 1:5
               allPts = [trialStruct([trialStruct.posture]==posture).PC]; 
               avg = mean(allPts);
               plot(1400,avg,'.','MarkerSize',30,'Color',pcmap(posture,:));
            end
            xticks([0 600 1200])
        end
        xlabel('trial'); ylabel('Posture Dim 1 (a.u.)')
        %title(dataset)
        set(gca, 'TickDir', 'out')
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_postureAxisVTrial.svg']));
        end
            
            
    end
    
