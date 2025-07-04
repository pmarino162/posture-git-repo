clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = "C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\rev 3 m 2 variance analysis";
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
    regularize = true;

%% Set up result struct
    resultStruct = struct('monkey','','dataset','','dPCATotalVar',[],'dPCAAbsVar',[],'dPCAPctVar',[],'PCAAbsVar',[],'PCAPctVar',[],'PCAPctVarOrderedLtoS',[]);
    resultStructInd = 1;
    
%% Main Loop  
    for datasetList = {'E20200318','N20171215','R20201020','R20201021'}%{'E20200318','E20210901','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221','E20210708'}
        %% Load Data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);

        %% Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams);
    
        %% Only use N00 and A90 for Nigel to match Rocky's postures
        if strcmpi(dataset,'N20171215')
            trajStruct = trajStruct([trajStruct.posture] ~= 2);
        end
        
        %% For earl reaching, keep postures 2 and 7, but relabel as 1 and 2 for color purposes
        if strcmpi(dataset,'E20210708')
            trajStruct = trajStruct(ismember([trajStruct.posture],[2,7]));
            for i = 1:size(trajStruct,2)
               if trajStruct(i).posture == 2
                   trajStruct(i).posture = 1;
               elseif trajStruct(i).posture == 7
                   trajStruct(i).posture = 2;
               end
            end
        end
              
        %% Get traj struct dimensions
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct); 
        
        %% Setup colormap (based on number of postures)
        [pcmap,tcmap,rainbow] = getColorMaps(numPostures);  

        %% Run dPCA
        [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADims(trajStruct,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps);   
                     

        %% Fit PCs
        [allTraj,mu] = collectAllAvgTraj(trajStruct);
        [PCs,~,~,~,explained,~] = pca(allTraj);
        

        %% Fill result struct
        numNeurons = size(allTraj,2);
        resultStruct(resultStructInd).monkey = dataset(1);
        resultStruct(resultStructInd).dataset = dataset;
        resultStruct(resultStructInd).totalVar = trace(cov(allTraj));
        resultStruct(resultStructInd).dPCAPctVar = [additionalVarExpl(postureDims(1)),additionalVarExpl(targetDims(1)),additionalVarExpl(targetDims(2))];
        resultStruct(resultStructInd).dPCAAbsVar = (resultStruct(resultStructInd).dPCAPctVar / 100) * resultStruct(resultStructInd).totalVar;
        resultStruct(resultStructInd).dPCANormAbsVar = ((resultStruct(resultStructInd).dPCAPctVar / 100) * resultStruct(resultStructInd).totalVar) / numNeurons;
        resultStruct(resultStructInd).normTotalVar = resultStruct(resultStructInd).totalVar / numNeurons;
        resultStruct(resultStructInd).eigenspec = explained;      
        
        
        %% Plot and save
        
        % Plot 1 - Absolute Variance
        figure;
        y = [resultStruct(resultStructInd).dPCAAbsVar, resultStruct(resultStructInd).totalVar];
        bar(y);
        ylabel("Absolute Variance")
        xticklabels({"Posture Dim 1", "Goal Dim 1", "Goal Dim 2", "Total"});
        title(['Monkey ',dataset(1)]);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_absVar.png']));
        end
        
        % Plot 2 - Pct Variance
        figure;
        y = resultStruct(resultStructInd).dPCAPctVar;
        bar(y);
        ylabel("Pct Variance")
        xticklabels({"Posture Dim 1", "Goal Dim 1", "Goal Dim 2"});
        title(['Monkey ',dataset(1)]);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_pctVar.png']));
        end
        
        % Plot 3 - Normalized Absolute Variance
        figure;
        y = [resultStruct(resultStructInd).dPCANormAbsVar, resultStruct(resultStructInd).normTotalVar];
        bar(y);
        ylabel("Absolute Variance Per Neuron")
        xticklabels({"Posture Dim 1", "Goal Dim 1", "Goal Dim 2", "Total"});
        title(['Monkey ',dataset(1)]);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_normAbsVar.png']));
        end
        
        % Plot 4 - Eigenspectrum
        figure;
        y = resultStruct(resultStructInd).eigenspec(1:20);
        bar(y);
        ylabel("Explained Variance by PC (%)")
        xlabel("PC");
        title(['Monkey ',dataset(1)]);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_eigenSpec.png']));
        end
        
        
        
        
        %% Update Result Struct Ind
        resultStructInd = resultStructInd + 1;
    end

    