clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

%     % BCI
%     [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForTrajDist(Data,'BCI');
%     trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
%     condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Step 2','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

    %Across posture BCI
    [Data] = loadEarlData20210901_20211210;
    [Data] = subselectForTrajDist(Data,'Multi Joint BCI');
    
    
%       %Iso force
%     [Data] = loadEarlData20200116_20211210();
%     [Data] = subselectForTrajDist(Data,'iso');
%     trialInclStates(1).trialName = {'IsometricForce_1D'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     
% %     Reaching
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForTrajDist(Data,'reaching');
%     trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%     %     Planning
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForTrajDist(Data,'planning');
%      trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Delay','first',250}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
 

%% Get timestamps dist and min number timestamps
    numTimestamps = [];
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allSmoothFR,2);
       numCondTraj = [numCondTraj,numTraj];
       for j = 1:numTraj
          numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
       end
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    
    figure
    histogram(numCondTraj)
    xlabel('Number of trials')
    ylabel('Number of conditions')
    
% Get minimum number of trials and timestamps
    [minNumTimestamps,i] = min(numTimestamps);
    [minNumCondTraj,i] = min(numCondTraj);


%% Do PCA on condition averages for visualization 
    avgSmoothFR = [trajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);       
    % Add Projections to trajStruct 
    for i = 1:size(trajStruct,2)
       trajStruct(i).PCA = trajStruct(i).avgSmoothFR.traj*coeff;
    end    
    
    

%% Do the heavy lifting
    % Create distance distributions
    resultStruct = struct('cond1P',[],'cond1T',[],'cond2P',[],'cond2T',[],'dist',[]);
    numCond = size(trajStruct,2);
    numDraws = 100; %Number of random draws from each condition
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    numPts = minNumTimestamps; %Number of points from each trajectory
    
    structInd = 1;
    for cond1 = 1:numCond
        cond1
        for cond2 = cond1:numCond
            resultStruct(structInd).cond1P = trajStruct(cond1).posture;
            resultStruct(structInd).cond1T = trajStruct(cond1).target;
            resultStruct(structInd).cond2P = trajStruct(cond2).posture;
            resultStruct(structInd).cond2T = trajStruct(cond2).target;
            for i = 1:numDraws
                %Create traj1
                numTraj1 = size(trajStruct(cond1).allSmoothFR,2);
                sampInd1 = randsample(numTraj1,numSample);
                traj1Struct = trajStruct(cond1).allSmoothFR(sampInd1);
                traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                %Create traj2
                if cond1 == cond2
                    numTrajRemaining = numTraj1 - numSample;
                    sampInd2 = randsample(numTrajRemaining,numSample);
                    remainingInd = setdiff(1:numTraj1,sampInd1);
                    sampInd2 = remainingInd(sampInd2);
                else
                    numTraj2 = size(trajStruct(cond2).allSmoothFR,2);
                    sampInd2 = randsample(numTraj2,numSample);
                end
                traj2Struct = trajStruct(cond2).allSmoothFR(sampInd2);
                traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                %Get alpha
                alpha = getOptimalAlpha(traj1,traj2,numPts);
                %alpha = traj1(1,:)-traj2(1,:);
                %Shift traj2 by alpha
                unshiftTraj2 = traj2;
                traj2 = traj2 + alpha;
                %Get residual distance
                dist = getMeanDist(traj1,traj2,numPts);
                resultStruct(structInd).dist(i) = dist;
                %Visualize trajectories
%                 visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist)
            end
            structInd = structInd + 1;
        end
    end
    
%% Collect lists of condition variables
cond1P = [resultStruct.cond1P]; cond1T = [resultStruct.cond1T]; cond2P = [resultStruct.cond2P]; cond2T = [resultStruct.cond2T];

%% Plot results
figure
cond1Vec = [3 3 3 3];
cond1Str = ['P',num2str(cond1Vec(1)),'T',num2str(cond1Vec(2)),' vs ','P',num2str(cond1Vec(3)),'T',num2str(cond1Vec(4))];
dist = resultStruct(cond1P==cond1Vec(1) & cond1T==cond1Vec(2) & cond2P==cond1Vec(3) & cond2T==cond1Vec(4)).dist;
histogram(dist)

hold on

cond2Vec = [1 3 3 3];
cond2Str = ['P',num2str(cond2Vec(1)),'T',num2str(cond2Vec(2)),' vs ','P',num2str(cond2Vec(3)),'T',num2str(cond2Vec(4))];
dist = resultStruct(cond1P==cond2Vec(1) & cond1T==cond2Vec(2) & cond2P==cond2Vec(3) & cond2T==cond2Vec(4)).dist;
histogram(dist)

legend(cond1Str,cond2Str)
xlabel('Mean Dist (Hz)')
ylabel('Count')

%% Check to see whether same condition distributions are overlapping 
for posture = 1:7
    for target = 1:8
        condVec = [posture target posture target];
        if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
            dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
            histogram(dist)
            hold on
        end
    end
end

xlabel('Mean Dist (Hz)')
ylabel('Count')

%% They are, so create chance dist out of them
chanceDist = [];
for posture = 1:7
    for target = 1:8
        condVec = [posture target posture target];
        if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
            dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
            chanceDist = [dist,chanceDist];
        end
    end
end

%% Compare conditions to chance dist

figure
histogram(chanceDist)
hold on

for target = 1:8
    for posture1 = 1:7
        for posture2 = posture1:7
            condVec = [posture1 target posture2 target];
            if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                histogram(dist)
            end
        end
    end
end

%% Combine all test dist
testDist = [];
for target = 1:8
    for posture1 = 1:7
        for posture2 = posture1+1:7
            condVec = [posture1 target posture2 target];
            if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                testDist = [dist,testDist];
            end
        end
    end
end

%% Combine all non test dist
nonTestDist = [];
storeCond = NaN(100000,4);
i = 1;
for posture1 = 1:7
    for target1 = 1:8
        for posture2 = posture1:7
            for target2 = target1+1:8
                condVec = [posture1 target1 posture2 target2];
                storeCond(i,:) = condVec;
                if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                    dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                    nonTestDist = [dist,nonTestDist];
                    i = i+1;
                end
                
            end
        end
    end
end

%% Lowest upper bound 
lowUB = [];
storeCond = NaN(100000,4);
i = 1;
for posture1 = 1:7
    for target1 = 1:8
        for posture2 = posture1
            for target2 = target1+1
                condVec = [posture1 target1 posture2 target2];
                if any(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4))
                    storeCond(i,:) = condVec;
                    dist = resultStruct(cond1P==condVec(1) & cond1T==condVec(2) & cond2P==condVec(3) & cond2T==condVec(4)).dist;
                    lowUB = [dist,lowUB];
                    i = i+1;
                end
                
            end
        end
    end
end



%% Plot chance and test
figure
histogram(chanceDist,'Normalization','probability')
hold on
histogram(testDist,'Normalization','probability')
histogram(nonTestDist,'Normalization','probability')
% histogram(lowUB,'Normalization','probability')
ax = gca;
minY = ax.YLim(1);
maxY = ax.YLim(2);
chanceMean = mean(chanceDist);
testMean = mean(testDist);
nonTestMean = mean(nonTestDist);
lowUBMean = mean(lowUB);
line([chanceMean chanceMean],[minY maxY],'Color','k','LineWidth',1.5)
line([testMean testMean],[minY maxY],'Color','k','LineWidth',1.5)
line([nonTestMean nonTestMean],[minY maxY],'Color','k','LineWidth',1.5)
% line([lowUBMean lowUBMean],[minY maxY],'Color','k','LineWidth',1.5)
xlabel('Mean Dist (Hz)')
ylabel('Probability')
legend('Within Condition Comparisons',['Comparisons w Different Postures',newline,'and Same Target'],'Comparisons w Different Targets') 
% legend('Within Condition Comparisons',['Comparisons w Different Postures',newline,'and Same Target'],'Comparisons w Different Targets','Comparisons with Adjacent Targets') 


% %% Test local functions
%     close all
%      
%     traj1 = [1 1;2 2];
% %     traj2 = traj1 + [1 1];
%     traj2 = [1 2; 2 4];
%     numPts = 2;
%     
%     figure
%     plot(traj1(:,1),traj1(:,2))
%     hold on
%     plot(traj2(:,1),traj2(:,2))
%     
%     getMeanDist(traj1,traj2,numPts)
%     alpha = getOptimalAlpha(traj1,traj2,numPts)
%     
%     figure
%     traj2 = traj2 + alpha;
%     plot(traj1(:,1),traj1(:,2))
%     hold on
%     plot(traj2(:,1),traj2(:,2))
%     
%     dist = getMeanDist(traj1,traj2,numPts)
%% Define local functions
    %Get optimal alpha
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end
    
    %Get mean dist
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end
    
    %Visualize traj
    function [] = visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist)
        traj1Proj = traj1*coeff;
        traj2Proj = traj2*coeff;
        unshiftTraj2Proj = unshiftTraj2*coeff;
        traj1AvgProj = trajStruct(cond1).PCA;
        traj2AvgProj = trajStruct(cond2).PCA;
        
        figure
        hold on
        %Traj 1
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'Color','r');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            plot3(traj1Proj(pt,1),traj1Proj(pt,2),traj1Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
        end
        %Unshifted Traj 2
        plot3(unshiftTraj2Proj(:,1),unshiftTraj2Proj(:,2),unshiftTraj2Proj(:,3),'Color','b');
        %Shifted Traj 2
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'Color','g');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            plot3(traj2Proj(pt,1),traj2Proj(pt,2),traj2Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
        end
        %Avg Traj 1
        plot3(traj1AvgProj(:,1),traj1AvgProj(:,2),traj1AvgProj(:,3),'Color','r','LineWidth',2);
        plot3(traj1AvgProj(1,1),traj1AvgProj(1,2),traj1AvgProj(1,3),'.','MarkerSize',10,'Color','r');
        %Avg Traj 2
        plot3(traj2AvgProj(:,1),traj2AvgProj(:,2),traj2AvgProj(:,3),'Color','b','LineWidth',2);
        plot3(traj2AvgProj(1,1),traj2AvgProj(1,2),traj2AvgProj(1,3),'.','MarkerSize',10,'Color','b');

        %Get condition info 
        cond1P = trajStruct(cond1).posture;
        cond1T = trajStruct(cond1).target;
        cond2P = trajStruct(cond2).posture;
        cond2T = trajStruct(cond2).target;
        
        
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(dist),'Hz'])
        close 
    end