clear; clc;

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load data
    dataset = 'N20190226';
    [Data,zScoreParams] = loadData(dataset);

%% Get Traj
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR','marker'};
    trialInclStates = struct('trialName','','inclStates',[]);
    trialInclStates(1).trialName = {'Nigel Dissociation'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    
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
    
%% Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);    
    
%% Plot Positions
    for posture = 1:2
        for target = 1:8
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgMarker.traj;
            if posture == 1
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
            else
                plot(traj(:,1),traj(:,2),'--','Color',tcmap(target,:));
            end
            hold on
        end
    end
    
%% Plot PSTH's
    figure
    for posture = 1:2
       for target = [1,5]
          traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.traj;
          time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.timestamps;
          for neuron = 1:55
              subplot(8,8,neuron)
              if posture == 1
                  plot(time,traj(:,neuron),'-','Color',tcmap(target,:));
              else
                  plot(time,traj(:,neuron),'--','Color',tcmap(target,:));
              end
              hold on
          end
       end
    end
    
%% Perform PCA
    X = [];
    for i = 1:size(trajStruct,2)
        X = vertcat(X,trajStruct(i).avgZSmoothFR.traj);
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(X);
    
%% Plot PC's
    figure
    PCs = [2,5,6];
    for posture = 1:2
        for target = [4,8]
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.traj;
            traj = traj*coeff;
            plot3(traj(:,PCs(1)),traj(:,PCs(2)),traj(:,PCs(3)),'Color',tcmap(target,:));
            hold on;
            plot3(traj(1,PCs(1)),traj(1,PCs(2)),traj(1,PCs(3)),'.','MarkerSize',20,'Color',tcmap(target,:));
        end
    end
    
    
    figure 
    for posture = 1:2
        for target = [1:8]
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.traj;
            traj = traj*coeff;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.timestamps;
            for PC = 1:10
              subplot(2,5,PC)
              if posture == 1
                  plot(time,traj(:,PC),'-','Color',tcmap(target,:));
              else
                  plot(time,traj(:,PC),'--','Color',tcmap(target,:));
              end
              hold on
            end
        end
    end