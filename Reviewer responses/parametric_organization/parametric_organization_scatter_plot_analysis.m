clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = "C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\Parametric organization";
    set(0, 'DefaultFigureRenderer', 'painters');

%% Set parameters
    numPCsToKeep = 10;     %Num PCs to project data into before analysis

%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(5);

%% Datasets to include in analysis
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319'};

%% Perform all measurements for all sessions. Save to resultStruct
    resultStruct = struct('dataset',[],'posturePoints',[],'posturePointsInPCs',[],'distances',[],'distancesInPCs',[]);
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
        % Keep only postures with all targets; get trajStructDims
        [postureList,~,targetList,~,~,~] = getTrajStructDimensions(trajStruct);
        [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList);
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        % Get posture signal, then posture points
        [minNumTimestamps] = getMinNumTimestamps(trajStruct);
        [pSig,tSig] = getPandTsig(trajStruct,'avgZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets);
        posturePoints = squeeze(mean(pSig,1,'omitnan'));
        % Project posture points onto PCs
        [PCs,~,~,~,explained,~] = pca(posturePoints);
        posturePointsInPCs = posturePoints*PCs;
        
        % Compute distances
        distances = struct('pDist',[],'nDist',[]);
        distancesStructInd = 1;
        for pDist = 1:4
            nDist = [];
            for startP = 1:4
                endP = startP + pDist;
                if endP <= 5
                    distance = norm(posturePoints(startP,:) - posturePoints(endP,:));
                    nDist = [nDist, distance];
                end
            end
            distances(distancesStructInd).pDist = pDist;
            distances(distancesStructInd).nDist = nDist;
            distancesStructInd = distancesStructInd + 1;
        end

        % Add to result
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).posturePoints = posturePoints;
        resultStruct(structInd).posturePointsInPCs = posturePointsInPCs;
        resultStruct(structInd).distances = distances;
        structInd = structInd + 1;
    end
    
%% Combine and plot
figure;
hold on;
degreesBtPostures = 15;
cmap = winter(4);
for pDist = 1:4
    allPDistPoints = [];
    for session = 1:length(bciDatasetList)
        allPDistPoints = [allPDistPoints, resultStruct(session).distances(pDist).nDist];
    end
    
    scatter(pDist*degreesBtPostures*ones(1,length(allPDistPoints)),allPDistPoints,60,cmap(pDist,:),"filled");%pcmap(pDist,:),"filled");
end
xlim([0,75])
set(gca,'FontSize',20, 'FontName', 'Arial')
xlabel('Change in arm angle (deg)')
ylabel('Neural Distance (a.u.)')
if saveFig
    saveas(gcf,fullfile(saveDir,'parametric_org.svg'));
end


%% Plot with box and whiskers
figure; 
hold on;
degreesBtPostures = 15;
positions = (1:4)*degreesBtPostures;  % x-positions for each pDist
cmap = winter(4);

% Preallocate arrays for boxplot
allData  = [];
groupIdx = [];

% Also store data in a cell array so we can scatter it later
pDistDataCell = cell(1,4);

% Gather data
for pDist = 1:4
    allPDistPoints = [];
    for session = 1:length(bciDatasetList)
        allPDistPoints = [allPDistPoints, resultStruct(session).distances(pDist).nDist];
    end
    
    % Save in a cell array for later plotting
    pDistDataCell{pDist} = allPDistPoints;
    
    % Also accumulate for boxplot
    allData  = [allData,  allPDistPoints];
    groupIdx = [groupIdx, pDist*ones(1,length(allPDistPoints))];
end

% Create the boxplot (default is median, box, whiskers)
hBox = boxplot(allData, groupIdx, ...
    'Positions', positions, 'Widths', 4, ...
    'Symbol','');
set(hBox, 'LineWidth', 2);

% Make the whiskers solid
set(findobj(gca,'Tag','Upper Whisker'),'LineStyle','-');
set(findobj(gca,'Tag','Lower Whisker'),'LineStyle','-');


% Plot the individual points with partial transparency
for pDist = 1:4
    scatter(positions(pDist)*ones(1,length(pDistDataCell{pDist})), ...
            pDistDataCell{pDist}, ...
            20, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.4);  % reduce alpha to de-emphasize points
end

% Axis and tick labels
set(gca, 'XTick', positions, ...
         'XTickLabel', {'15','30','45','60'}, ...
         'FontSize', 20, 'FontName', 'Arial');
xlim([0, 75]);
xlabel('Change in arm angle (deg)');
ylabel('Neural Distance (a.u.)');
if saveFig
    saveas(gcf,fullfile(saveDir,'parametric_org_box_and_whisker.svg'));
end