clear; clc; clf; close all; 

%% Load Data 
    %Earl - 3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);
    numPostures = 5;
    numTargets = 8;

%% Create trajStruct
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'allChannelSmoothedFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
        trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
        trialInclStates(1).addTimeToBeginning = {-100};
        trialInclStates(1).addTimeToEnd = {0};
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates);
    
%% Use all data to perform posture, target, and PxT LDA
    %Preallocate
    numDims = 87;
    totalLength = 0;
    for i = 1:size(trajStruct,2)
       numTraj = size(trajStruct(i).allAllChannelSmoothedFR,2);
       for j = 1:numTraj
           numSteps = size(trajStruct(i).allAllChannelSmoothedFR(j).traj,1);
           totalLength = totalLength + numSteps;
       end 
    end
    allTraj = zeros(totalLength,numDims); 
    allPostureLabels = zeros(totalLength,1); allTargetLabels = zeros(totalLength,1); allPxTLabels = zeros(totalLength,1);
    
    %Fill 
    curRow = 1;
    for i = 1:size(trajStruct,2)
        i
       target = trajStruct(i).target;
       posture = trajStruct(i).posture;
       numTraj = size(trajStruct(i).allAllChannelSmoothedFR,2);
       for j = 1:numTraj
            %Trajectory
            traj = trajStruct(i).allAllChannelSmoothedFR(j).traj;
            numSteps = size(traj,1);
            allTraj(curRow:curRow+numSteps-1,:) = traj;
            %Labels
            allTargetLabels(curRow:curRow+numSteps-1,:) = ones(numSteps,1).*target;
            allPostureLabels(curRow:curRow+numSteps-1,:) = ones(numSteps,1).*posture;
            allPxTLabels(curRow:curRow+numSteps-1,:) = ones(numSteps,1).*str2double([num2str(posture),num2str(target)]);
            curRow = curRow + numSteps;
        end
    end
    postureLDA = fisherLDA(allTraj, allPostureLabels);
    [postureLDA,~] = qr(postureLDA);
    postureLDA = postureLDA(:,1:4);
    
    targetLDA = fisherLDA(allTraj, allTargetLabels);
    [targetLDA,~] = qr(targetLDA);
    targetLDA = targetLDA(:,1:7);
    
    [postTargOrth,~] = qr([postureLDA(:,1),targetLDA(:,1:2)]);
    postTargOrth = postTargOrth(:,1:3);
    
    nullSpace = null(postTargOrth');
    nullProj = allTraj*nullSpace;
    nullSpacePC = pca(nullProj);
    nullSpace = nullSpace*nullSpacePC;
        
    pxtLDA = fisherLDA(allTraj, allPxTLabels);
    [pxtLDA,~] = qr(pxtLDA);
    pxtLDA = pxtLDA(:,1:numDims-1);
    
%% Get Neutral Decoder Workspace
    N00CRinv = N00DecoderParams.CRinv; N00W = N00DecoderParams.W;
    N00C = N00GPFAParams.estParams.C; N00d = N00DecoderParams.d;
    [N00U, N00D, N00V] = svd(N00C, 0); 
    N00Decoder = N00W*N00V*inv(N00D)*N00U';
    N00Decoder = N00Decoder';
    [N00DecodePlane,~] = qr(N00Decoder);
    N00DecodePlane = N00DecodePlane(:,1:2);
    
    E30CRinv = E30DecoderParams.CRinv; E30W = E30DecoderParams.W;
    E30C = E30GPFAParams.estParams.C; E30d = E30DecoderParams.d;
    [E30U, E30D, E30V] = svd(E30C, 0); 
    E30Decoder = E30W*E30V*inv(E30D)*E30U';
    E30Decoder = E30Decoder';
    [E30DecodePlane,~] = qr(E30Decoder);
    E30DecodePlane = E30DecodePlane(:,1:2);
    
    I30CRinv = I30DecoderParams.CRinv; I30W = I30DecoderParams.W;
    I30C = I30GPFAParams.estParams.C; I30d = I30DecoderParams.d;
    [I30U, I30D, I30V] = svd(I30C, 0); 
    I30Decoder = I30W*I30V*inv(I30D)*I30U';
    I30Decoder = I30Decoder';
    [I30DecodePlane,~] = qr(I30Decoder);
    I30DecodePlane = I30DecodePlane(:,1:2);
    
%% Add Projections to Traj Struct
    for i = 1:size(trajStruct,2)
       trajStruct(i).postureLDA = trajStruct(i).avgAllChannelSmoothedFR.traj*postureLDA;
       trajStruct(i).targetLDA = trajStruct(i).avgAllChannelSmoothedFR.traj*targetLDA;
       trajStruct(i).postTargOrth = trajStruct(i).avgAllChannelSmoothedFR.traj*postTargOrth;
       trajStruct(i).PxTLDA = trajStruct(i).avgAllChannelSmoothedFR.traj*pxtLDA;
       trajStruct(i).null = trajStruct(i).avgAllChannelSmoothedFR.traj*nullSpace;
       trajStruct(i).N00Decoder = trajStruct(i).avgAllChannelSmoothedFR.traj.*(45/1000)*N00Decoder;
       trajStruct(i).I30Decoder = trajStruct(i).avgAllChannelSmoothedFR.traj.*(45/1000)*I30Decoder;
       trajStruct(i).E30Decoder = trajStruct(i).avgAllChannelSmoothedFR.traj.*(45/1000)*E30Decoder;
    end
    
%% RSV
%Set up spaceIDs
spaceIDs = struct('ID',[],'spaceName',[],'dirs',[]);
spaceIDs(1).ID = 1; spaceIDs(1).spaceName = 'N00DecodePlane'; spaceIDs(1).dirs = N00DecodePlane;
spaceIDs(2).ID = 2; spaceIDs(2).spaceName = 'I30DecodePlane'; spaceIDs(2).dirs = I30DecodePlane;
spaceIDs(3).ID = 3; spaceIDs(3).spaceName = 'E30DecodePlane'; spaceIDs(3).dirs = E30DecodePlane;
spaceIDs(4).ID = 4; spaceIDs(4).spaceName = 'postureLDA1'; spaceIDs(4).dirs = postureLDA(:,1);
spaceIDs(5).ID = 5; spaceIDs(5).spaceName = 'postureLDA2'; spaceIDs(5).dirs = postureLDA(:,2);
spaceIDs(6).ID = 6; spaceIDs(6).spaceName = 'targetLDA'; spaceIDs(6).dirs = targetLDA(:,1:2);

%Compute RSV and Principal Angles
RSV = struct('space1','','space2','','RSV',[],'prinAngles',[],'vectorLengths',[]);
structInd = 1;
numSpaces = size(spaceIDs,2);
heatmapMat = zeros(numSpaces);
for space1ID = 1:numSpaces
    for space2ID = 1:numSpaces
        space1Dirs = spaceIDs(space1ID).dirs;
        space2Dirs = spaceIDs(space2ID).dirs;
        RSV(structInd).space1 = spaceIDs(space1ID).spaceName;
        RSV(structInd).space2 = spaceIDs(space2ID).spaceName;
        RSV(structInd).RSV = sum(var(allTraj*space1Dirs*space1Dirs'*space2Dirs))/sum(var(allTraj*space1Dirs));
        RSV(structInd).prinAngles = rad2deg(prinangle(space1Dirs,space2Dirs));
        vectorLengths = zeros(1,size(space1Dirs,2));
        for i=1:size(vectorLengths,2)
            vectorLengths(i) = norm(space1Dirs(:,i)'*space2Dirs);
        end
        RSV(structInd).vectorLengths = vectorLengths;
        heatmapMat(space1ID,space2ID) = sum(var(allTraj*space1Dirs*space1Dirs'*space2Dirs))/sum(var(allTraj*space1Dirs));
        structInd = structInd + 1;
    end
end

%% Plot table 
hold off
xvalues = {spaceIDs.spaceName};
yvalues = {spaceIDs.spaceName};
h = heatmap(xvalues,yvalues,heatmapMat);

%% Length of posture axes in decode plane 
axis1 = postureLDA(:,1); axis2 = N00Decoder(:,1);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)
projLength = dot(axis1,axis2)
postureLDA(:,1)'*(45/1000)*N00Decoder;

axis1 = postureLDA(:,2); axis2 = N00Decoder(:,1);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)
projLength = dot(axis1,axis2)





%% Get Angles and Projection Lengths
projLengths = struct('proj','','length',[]);
%Posture LDA1 onto N00 Decode Plane

%Posture LDA2 onto N00 Decode Plane 



axis1 = postureLDA(:,1); axis2 = targetLDA(:,1);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)

axis1 = postureLDA(:,1); axis2 = targetLDA(:,2);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)

axis1 = postureLDA(:,1); axis2 = N00Decoder(:,1);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)

axis1 = targetLDA(:,1); axis2 = N00Decoder(:,1);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)

axis1 = postureLDA(:,2); axis2 = N00Decoder(:,2);
theta = acosd(dot(axis1,axis2)/(norm(axis1)*norm(axis2)))
projLength = dot(axis1,axis2)/norm(axis2)

%% Plot projection lengths
close all
figure
RSVInds = [19,20,21,24,25,26,27,30];
lengths = [RSV(RSVInds).vectorLengths];
bar(lengths);
ylabel('Length')
space1IDs = {RSV(RSVInds).space1};
space2IDs = {RSV(RSVInds).space2};
labels = cell(1,length(RSVInds));
for i=1:length(RSVInds)
   labels{i} = [space1IDs{i},' onto ',space2IDs{i}]; 
    
end
xticklabels(labels)
xtickangle(45)
% hold off
% xvalues = {spaceIDs.spaceName};
% yvalues = {spaceIDs.spaceName};
% h = heatmap(xvalues,yvalues,heatmapMat);

%% Plot
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
%N00 Decoder    
figure 
xDim = 1; yDim = 2;
for posture = [1,3,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).N00Decoder;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('N00 Decoder')
%E30 Decoder
figure 
xDim = 1; yDim = 2;
for posture = [1,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).E30Decoder;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('E30 Decoder')
%I30 Decoder
figure 
xDim = 1; yDim = 2;
for posture = [1,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).I30Decoder;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('I30 Decoder')
%Target LDA
figure 
xDim = 1; yDim = 2;
for posture = [1,5]
    for target = 1:numTargets
        traj = trajStruct(find([trajStruct.posture]==posture & [trajStruct.target]==target)).targetLDA;
        plot(traj(:,xDim),traj(:,yDim),'Color',cmap(posture,:),'LineWidth',2); 
        hold on
    end
end
title('Target LDA')