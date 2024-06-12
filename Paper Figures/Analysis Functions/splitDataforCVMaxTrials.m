function [trajStruct1,trajStruct2] = splitDataforCVMaxTrials(trajStruct,binWidth)       
    % Alternative method for splitting data that uses the max possible
    % number of trials for each condition (used for assessing
    % shoulder/elbow space alignment)

    %Split trajStruct into two groups for cross-validation
    trajStruct1 = trajStruct;
    trajStruct2 = trajStruct;
    for j = 1:size(trajStruct,2)
        numTraj = size(trajStruct(j).allZSmoothFR,2);
        numSample = floor(numTraj/2); %Number of trajectories in each draw
        sampInd1 = randsample(numTraj,numSample);     
        trajStruct1(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd1);
        [trajStruct1(j).avgZSmoothFR.traj,trajStruct1(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct1(j),'zSmoothFR',binWidth);                       
        trajStruct1(j).allPCA = trajStruct(j).allPCA(sampInd1);
        [trajStruct1(j).avgPCA.traj,trajStruct1(j).avgPCA.timestamps,trajStruct1(j).avgPCA.CI95] = getAvgTraj(trajStruct1(j),'PCA',binWidth); 
        
        numTrajRemaining = numTraj - numSample;
        sampInd2 = randsample(numTrajRemaining,numSample);
        remainingInd = setdiff(1:numTraj,sampInd1);
        sampInd2 = remainingInd(sampInd2);
        trajStruct2(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd2);
        [trajStruct2(j).avgZSmoothFR.traj,trajStruct2(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct2(j),'zSmoothFR',binWidth);  
        trajStruct2(j).allPCA = trajStruct(j).allPCA(sampInd2);
        [trajStruct2(j).avgPCA.traj,trajStruct2(j).avgPCA.timestamps,trajStruct2(j).avgPCA.CI95] = getAvgTraj(trajStruct2(j),'PCA',binWidth); 
                
    end   
end