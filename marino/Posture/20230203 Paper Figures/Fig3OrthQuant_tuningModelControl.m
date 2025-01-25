clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20221123\Figures\Figure 2\Alignment Indices\20230125\fullSpace';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    scmap = customRainbow;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'vPP',[],'vPT',[],'vTT',[],'vTP',[],'APT',[],'ATP',[],'nullAPT',[],'nullATP',[]);
    structInd = 1;

%% Generate Neural and Kinematic Data
    posLag = 125;
    velLag = 75;
    [kinStruct,neuronStruct,Data] = generateNeuralAndKinData_20230210(posLag,velLag);

%% Create whole trial, reaching, and holding traj structs
    dataset = 'E20230210';
    %Parameters
    trajFields = {'allChannelSmoothedFR','marker','markerVel'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'Center Out'};

    %Reaching
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    
                trialInclStates(1).trialName = {'Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Hold','first',0}}; 
    

    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false);
        
       
%% Run Loop for each dataset
       
        %Get posture, target, and overall spaces
            %Get min timesteps
            numTimestamps = [];
            for i = 1:size(trajStruct,2)
                numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
            end
            [minNumTimestamps,i] = min(numTimestamps);
            
            %Get posture & target lists
            postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
            targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
            numChannels = size(trajStruct(1).avgSmoothFR.traj,2); numConditions = size(trajStruct,2);
            
           
            
            %Update posture list
            postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
            
            %Form X
            X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
            postureInd = 1;
            for posture = postureList
                targetInd = 1;
                for target = targetList
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                    X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
                    targetInd = targetInd + 1;
                end
                postureInd = postureInd+1;
            end
            
            %Compute within space variances
            allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
            [allPCs,~,~,~,explained,allMu] = pca(allTraj);
            
            % Do marginalizations of X
            %Xdims: 1=time, 2=target, 3=posture, 4=channel
            %Condition-invariant
            CIMargOffset = mean(X,[1 2 3],'omitnan');
            CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
            %Posture and Target Traj
            targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
            postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
        
            %Compute PCs of marginalizations
            CIMargTraj = squeeze(CIMargTraj);
                CIMargOffset = squeeze(CIMargOffset);
                [CIPCs,~,~,~,explainedCI,CIMu] = pca(CIMargTraj); 
            targetMargTraj = squeeze(targetMargTraj);
                targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
                [targetPCs,~,~,~,explainedT,targetMargMu] = pca(targetMargTraj); 
            postureMargTraj = squeeze(postureMargTraj);
                postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
                [posturePCs,~,~,~,explainedP,postureMargMu] = pca(postureMargTraj); 
            
            %Define num posture dims and num target dims
            numPdims = 2; numTdims = 2;


        %Compute within and cross-variances
        CP = cov(postureMargTraj);
        DP = posturePCs(:,1:numPdims);
        CT = cov(targetMargTraj);
        DT = targetPCs(:,1:numTdims);
        
        
        if numPdims==2 && numTdims==2
            vPP = diag(DP'*CP*DP)./trace(CP);
            vPT = diag(DT'*CP*DT)./trace(CP);
            vTT = diag(DT'*CT*DT)./trace(CT);
            vTP = diag(DP'*CT*DP)./trace(CT);
            APT = trace(DT'*CP*DT)/trace(DP'*CP*DP);
            ATP = trace(DP'*CT*DP)/trace(DT'*CT*DT);
        elseif numPdims==2 && numTdims==1 %target var captured by "best" posture dim
            [bestPdims,~,~,~,explained_bestP,~] = pca(targetMargTraj*DP); 
            vPP = diag(DP'*CP*DP)./trace(CP);
            vPP = vPP(1);
            vPT = diag(DT'*CP*DT)./(vPP*trace(CP));
            vTT = diag(DT'*CT*DT)./trace(CT);
            vTP = var(targetMargTraj*DP*bestPdims(:,1))./trace(CT); 
            
            APT = diag(DT'*CP*DT)./(vPP*trace(CP));
            ATP =  var(targetMargTraj*DP*bestPdims(:,1))/trace(DT'*CT*DT);
        elseif numPdims==1 && numTdims==2%posture var captured by "best" target dim
            [bestTdims,~,~,~,explained_bestT,~] = pca(postureMargTraj*DT); 
            vTT = diag(DT'*CT*DT)./trace(CT);
            vTT = vTT(1);           
            vTP = diag(DP'*CT*DP)./(vTT*trace(CT));
            vPP = diag(DP'*CP*DP)./trace(CP);   
            vPT = var(postureMargTraj*DT*bestTdims(:,1))./trace(CP); 

            APT = var(postureMargTraj*DT*bestTdims(:,1))/trace(DP'*CP*DP);            
            ATP =  diag(DP'*CT*DP)./(vTT*trace(CT));   
        end
        resultStruct(structInd).vPP = vPP;
        resultStruct(structInd).vPT = vPT;
        resultStruct(structInd).vTT = vTT;
        resultStruct(structInd).vTP = vTP;
        resultStruct(structInd).APT = APT;
        resultStruct(structInd).ATP = ATP;

        %Compute chance distribution
        numReps = 1000;
        nullAPT = NaN(1,numReps);
        nullATP = NaN(1,numReps);
        for i = 1:numReps
            if numPdims == 2 && numTdims == 2
                v1 = normrnd(0,1,numChannels,2);
                v2 = normrnd(0,1,numChannels,2);
                DR = orth([v1,v2]);
                APT = trace(DR'*CP*DR)/trace(DP'*CP*DP);
                ATP = trace(DR'*CT*DR)/trace(DT'*CT*DT);
                nullAPT(i) = APT;
                nullATP(i) = ATP;
            else
                v = normrnd(0,1,numChannels,2);
                DR = orth(v);
                APT = trace(DR'*CP*DR)/trace(DP'*CP*DP);
                ATP = trace(DR'*CT*DR)/trace(DT'*CT*DT);
                nullAPT(i) = APT;
                nullATP(i) = ATP;
            end
        end

        %Get dot products between axes
        resultStruct(structInd).dT1P1 = dot(DT(:,1),DP(:,1));
        if numTdims > 1
            resultStruct(structInd).dT2P1 = dot(DT(:,2),DP(:,1));
            if numPdims > 1
                resultStruct(structInd).dT1P2 = dot(DT(:,1),DP(:,2));
                resultStruct(structInd).dT2P2 = dot(DT(:,2),DP(:,2));
            end
        end
        if numPdims > 1 && numTdims < 2
            resultStruct(structInd).dT1P2 = dot(DT(:,1),DP(:,2));
        end


            %Add to resultStruct
            resultStruct(structInd).animal = dataset(1);
            resultStruct(structInd).dataset = dataset;
            resultStruct(structInd).nullAPT = nullAPT;
            resultStruct(structInd).nullATP = nullATP;
            structInd = structInd + 1;

    
%% Plot Results
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).dataset(1);
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    fs = 12;
    f=figure; hold on;
    %f.Position = [200 200 400 600];
    monkeyInd = 1;
    xtickList = []; xtickListInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       APT = horzcat(tempResultStruct.APT);
       ATP = horzcat(tempResultStruct.ATP);
       nullAPT = horzcat(tempResultStruct.nullAPT);
       nullAPTmean = mean(nullAPT);
       nullAPTstd = std(nullAPT);
       nullATP = horzcat(tempResultStruct.nullATP);
       nullATPmean = mean(nullATP);
       nullATPstd = std(nullATP);
       numSamp = size(APT,2);
       %Plot Data and errorbars
       offset = (numMonkeys*2)/(5*numMonkeys-1);
       alpha = 0.3;
       %APT
           scatter(monkeyInd*2*ones(1,numSamp)-1.5*offset,APT,100,scmap(6,:),'filled','MarkerFaceAlpha',0.5);
           e = errorbar(monkeyInd*2-0.5*offset,nullAPTmean,nullAPTstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k');      
           set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
           [h,p,ci,stats] = ttest2(APT,nullAPT);
           text(monkeyInd*2-1*offset,0.9,['p=',num2str(p)]);
       %ATP
           scatter(monkeyInd*2*ones(1,numSamp)+0.5*offset,ATP,100,scmap(6,:),'filled','MarkerFaceAlpha',0.5);
           e = errorbar(monkeyInd*2+1.5*offset,nullATPmean,nullATPstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k');
           set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
           [h,p,ci,stats] = ttest2(ATP,nullATP);
           text(monkeyInd*2+1*offset,0.9,['p=',num2str(p)]);
       %Update xtickList
       xtickList(xtickListInd) = monkeyInd*2 - 1.5*offset;
       xtickList(xtickListInd+1) = monkeyInd*2 - 0.5*offset;
       xtickList(xtickListInd+2) = monkeyInd*2 + 0.5*offset;
       xtickList(xtickListInd+3) = monkeyInd*2 + 1.5*offset;
       xtickListInd = xtickListInd + 4;
       monkeyInd = monkeyInd + 1;
    end
    
    ax = gca;
    ax.YLim = [0 1]; ax.XLim = [1 numMonkeys*2+1];
    yticks([0 0.5 1]); 
   
    %xtickLabelNames = {['Post',newline,'in',newline,'Targ'],['Targ',newline,'in',newline,'Post'],'Random'};
    
    row1 = {'Post Sig','Rand','Targ Sig','Rand'};
    row2 = {'in','','in',''};
    row3 = {'Targ','','Post',''};
    row4 = {'Space','','Space',''};
    
    labelArray = [row1;row2;row3;row4];
    if numMonkeys ==2 
        labelArray = [row1 row1 ;row2 row2; row3 row3; row4 row4];
    elseif numMonkeys == 3
        labelArray = [row1 row1 row1; row2 row2 row2; row3 row3 row3; row4 row4 row4];
    end
    labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center
    tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\\newline%s\n', labelArray{:}));

    xticks(xtickList);
    xticklabels(tickLabels)
    set(gca, 'TickDir', 'out')
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('Normalized Variance Captured')
    %xlabel('Monkey')
    if saveFig
        saveas(gcf,fullfile(saveDir,[task,'_alignmentInd.svg']));
    end
    
    
    

