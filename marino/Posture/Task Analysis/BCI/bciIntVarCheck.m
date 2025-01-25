clear; clc; clf; close all;

%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Load Data   
    binWidth = 25;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);

%% Get Traj Struct   
    trajFields = {'decoderGPFATraj'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
            condFields = {{'posture','conditionData','postureID'},{'target','targetData','target1ID'}};
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            trialInclStates(1).inclStates = {'Step 1'};
            trialInclStates(1).inclOccurrence = {'last'};
                trialInclStates(1).addTimeToBeginning = {0};
                trialInclStates(1).addTimeToEnd = {0};   
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,'matchConditions',true);

%% Add workspace projections
    trajStruct = trajStruct([trajStruct.posture]==3);
    W = N00DecoderParams.W;
    W = orth(W')';
    for i = 1:size(trajStruct,2)
       traj = trajStruct(i).avgDecoderGPFATraj.traj;
       trajStruct(i).avgDecoderWorkTraj = (W*traj')';
    end
    nullBasis = null(W);

%% Project Condition averages into null space of posture axis
    %Vertically concatenate trial averages
    allAvgs = [];
    for i = 1:size(trajStruct,2)
        traj = trajStruct(i).avgDecoderGPFATraj.traj;
        allAvgs = vertcat(allAvgs,traj);
    end
    nullProj = allAvgs*nullBasis;
    [nullBasisPC,score,latent,tsquared,explained,mu] = pca(nullProj);
    nullBasis = nullBasis*nullBasisPC;

%% Add workspace projections
    for i = 1:size(trajStruct,2)
       traj = trajStruct(i).avgDecoderGPFATraj.traj;
       trajStruct(i).avgDecoderNullProj = traj*nullBasis;
    end
    
%% Plot neural activity 
    figure
    for target = 1:8
        traj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj;
        plot(traj(:,1),traj(:,2));
        hold on
        
    end

    figure
    
    for target = 1:8
        time = trajStruct([trajStruct.target]==target).avgDecoderGPFATraj.timestamps;
        targTraj = trajStruct([trajStruct.target]==target).avgDecoderWorkTraj;
        nullTraj = trajStruct([trajStruct.target]==target).avgDecoderNullProj;
        
        for dim = 1:2
            subplot(2,5,dim)
            plot(time,targTraj(:,dim))
            hold on  
        end
        for dim = 3:10
            subplot(2,5,dim)
            plot(time,nullTraj(:,dim-2))
            hold on  
        end
        
    end

    
    
    

    
%% Get Variance Explained in condition averages
    %Get total variance in condition averages
        totalVar = sum(var(allAvgs));

    %Get variance explained in each PCA dim
        workVarExpl = (var(allAvgs*W')/totalVar)*100;
        nullVarExpl = (var(allAvgs*nullBasis)/totalVar)*100;
        