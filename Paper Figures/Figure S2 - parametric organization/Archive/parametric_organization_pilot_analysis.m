clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Reviewer responses\Analyses\Parametric organization';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Set parameters
    numPCsToKeep = 10;     %Num PCs to project data into before analysis

%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(5);

%% Datasets to include in analysis
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319'};

%% Perform all measurements for all sessions. Save to resultStruct
    resultStruct = struct('dataset',[],'distances',[],'posturePoints',[],'posturePointsInPCs',[]);
    structInd = 1;
    
    for datasetList = bciDatasetList
        % Load data
        dataset = datasetList{1,1};
        monkey = dataset(1);
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        % Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParams(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams); 
        %Project all data down to top PCs
        %[allTraj] = collectAllAvgTraj(trajStruct);
        %[trajStruct] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep);
        % Keep only postures with all targets; get trajStructDims
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList);
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        % Get posture signal, then posture points
        [minNumTimestamps] = getMinNumTimestamps(trajStruct);
        [pSig,tSig] = getPandTsig(trajStruct,'avgZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets);
        posturePoints = squeeze(mean(pSig,1,'omitnan'));
        % Compute distances
        distances = zeros(1,numPostures-1);
        for i = 1:numPostures-1
            distances(i) = norm(posturePoints(i,:) - posturePoints(i+1,:));
        end
        % Project posture points onto PCs
        [PCs,~,~,~,explained,~] = pca(posturePoints);
        posturePointsInPCs = posturePoints*PCs;
        % Add to result
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).distances = distances;
        resultStruct(structInd).posturePoints = posturePoints;
        resultStruct(structInd).posturePointsInPCs = posturePointsInPCs;
        structInd = structInd + 1;
    end
    
%% Plot - per-session PCs and bars
numSessions = size(resultStruct,2);

f = figure; 
f.Position = [100, 100, 1000, 1000];
plotInd = 1;

% Precompute limits for the first column scatter plots
allPC1 = [];
allPC2 = [];
for session = 1:numSessions
    posturePointsInPCs = resultStruct(session).posturePointsInPCs;
    allPC1 = [allPC1; posturePointsInPCs(:,1)];
    allPC2 = [allPC2; posturePointsInPCs(:,2)];
end
xLimits = [min(allPC1), max(allPC1)] + [-0.1, 0.1] * (max(allPC1) - min(allPC1));
yLimits = [min(allPC2), max(allPC2)] + [-0.1, 0.1] * (max(allPC2) - min(allPC2));

for session = 1:numSessions
    posturePointsInPCs = resultStruct(session).posturePointsInPCs;
    distances = resultStruct(session).distances;

    % First column: scatter plot (PC1 vs PC2)
    subplot(numSessions, 3, plotInd)
    for posture = postureList
        plot(posturePointsInPCs(posture, 1), posturePointsInPCs(posture, 2), '.', ...
             'MarkerSize', 20, 'Color', pcmap(posture, :));
        hold on;
    end
    xlabel('Posture PC 1')
    ylabel('Posture PC 2')
    xlim(xLimits); % Apply shared axis limits
    ylim(yLimits);
    plotInd = plotInd + 1;

    % Second column: scatter plot (PC3 vs PC4)
    subplot(numSessions, 3, plotInd)
    for posture = postureList
        plot(posturePointsInPCs(posture, 3), posturePointsInPCs(posture, 4), '.', ...
             'MarkerSize', 20, 'Color', pcmap(posture, :));
        hold on;
    end
    xlabel('Posture PC 3')
    ylabel('Posture PC 4')
    xlim(xLimits); % Match x-axis limits to first column
    ylim(yLimits); % Match y-axis limits to first column
    plotInd = plotInd + 1;

    % Third column: bar plot
    subplot(numSessions, 3, plotInd)
    bar([1:4], distances);
    xticks([1:4])
    xticklabels({'1->2', '2->3', '3->4', '4->5'});
    xlabel('Posture Pair')
    ylabel('Eucl. dist.')
    plotInd = plotInd + 1;

end

if saveFig
    saveas(gcf,fullfile(saveDir,'individual_sessions.png'));
end


%% Plot - Overall summary bar plot
figure; hold on;

resultMat = transpose(reshape([resultStruct.distances],4,numSessions));
means = mean(resultMat,1);
stdDev = std(resultMat);
stdErr = stdDev./(sqrt(size(resultMat,1)));

bar([1:4],means);

er = errorbar([1:4]+.05,means,stdErr);
er.Color = [1,0,0];
er.LineStyle = 'none';
er.LineWidth = 1;

for comparison = 1:4
    plot(comparison-.05, resultMat(:,comparison),'.k','MarkerSize',15);
end

xticks([1:4])
xticklabels({'1->2', '2->3', '3->4', '4->5'});
xlabel('Posture Pair')
ylabel('Euclidean distance in space of top 10 PCs')

if saveFig
    saveas(gcf,fullfile(saveDir,'summary_bar.png'));
end