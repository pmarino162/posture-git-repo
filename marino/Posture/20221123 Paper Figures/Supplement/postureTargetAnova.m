clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Supplement\Individual Neurons Are Mixed';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    cmap = customRainbow([3,4,5,6],:);
        
%% Run loop for each dataset    
    for datasetList = {'N20190226','R20200221'}
        %{'E20200317','N20171215','R20201020','E20200116','E20210706','N20190226','R20200221'}%{'R20201020'}%{'E20200116' 'E20210706','N20180221'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            %Iso
            case 'E20200116'
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
            %Reaching
            case 'E20210706'
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        %Keep postures 1-4 only for earl
        switch dataset
            case {'E20210706'}
                    trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4,5,7]));
        end
        
        %Get 1 observation per trial
        for i = 1:size(trajStruct,2)
            %Individual trials
            for j = 1:size(trajStruct(i).allSmoothFR,2)
               trajStruct(i).allSmoothFR(j).trialAvg = mean(trajStruct(i).allSmoothFR(j).traj); 
            end
        end
        
        %Get minNumber of reps
        numTrials = [];
        for i = 1:size(trajStruct,2)
            numTrials(i) = size(trajStruct(i).allSmoothFR,2);
        end
        [minNumTrials,i] = min(numTrials);
               
        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);
                       
        %Format Data for Anova
        dataTable = NaN(minNumTrials*numPostures,numTargets,numChannels);
        postureInd = 1;
        tableInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                allSmoothFR = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR;
                trialAvg = vertcat(allSmoothFR.trialAvg);
                dataTable(tableInd:tableInd+minNumTrials-1,targetInd,:) = datasample(trialAvg,minNumTrials,1,'Replace',false);
                targetInd = targetInd + 1;
            end
            tableInd = tableInd + minNumTrials;
            postureInd = postureInd + 1;
        end
        
        %Perform Anova 
        resultStruct = struct('channel',[],'dataTable',[],'p',[],'tbl',[],'stats',[]);
        structInd = 1;
        for channel = 1:numChannels
            resultStruct(structInd).channel = channel;
            resultStruct(structInd).dataTable = dataTable(:,:,channel);
            [p,tbl,stats] = anova2(resultStruct(structInd).dataTable,minNumTrials,'off');
            resultStruct(structInd).p = p;
            resultStruct(structInd).tbl = tbl;
            resultStruct(structInd).stats = stats;
            structInd = structInd + 1;
        end
        
        %Plot Pie Chart
        numTarget = 0; numPosture = 0; numBoth = 0; numNeither = 0;
        for channel = 1:numChannels
            targetP = resultStruct(channel).p(1,1);
            postureP = resultStruct(channel).p(1,2);
            if postureP < 0.05 
                if targetP < 0.05
                    numBoth = numBoth + 1;
                else
                    numPosture = numPosture + 1;
                end
            end
            if postureP > 0.05 && targetP < 0.05 
                numTarget = numTarget + 1;
            end
            if postureP > 0.05 && targetP > 0.05
                numNeither = numNeither + 1;
            end
        end
        
        figure; fs = 14;
        
        pieValues = [numTarget,numPosture,numBoth,numNeither];
        txt = {'Target Only: ';'Posture Only: ';'Mixed: ';'Neither: '}; 
        neuronCounts = {[' (',num2str(numTarget),')']; [' (',num2str(numPosture),')'];...
            [' (',num2str(numBoth),')'];[' (',num2str(numNeither),')']};
        
        p = pie(pieValues(pieValues~=0));
        colormap(cmap(pieValues~=0,:));
        
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 

        combinedtxt = strcat(txt(pieValues~=0,1),percentValues,neuronCounts(pieValues~=0,1)); 
        
        for i = 1:size(pText,1)
            pText(i).String = combinedtxt(i);
        end
        
        title(dataset);
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_tuningPieChart.svg']));
        end
        
    end