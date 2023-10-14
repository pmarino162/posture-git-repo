function [trajStruct1,trajStruct2] = splitDataforCV(trajStruct,minNumCondTrials,binWidth)       

    %Split trajStruct into two groups for cross-validation
    numSample = floor(minNumCondTrials/2); %Number of trajectories in each draw
    trajStruct1 = trajStruct;
    trajStruct2 = trajStruct;
    for j = 1:size(trajStruct,2)
        numTraj = size(trajStruct(j).allZSmoothFR,2);
        sampInd1 = randsample(numTraj,numSample);
        trajStruct1(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd1);
        [trajStruct1(j).avgZSmoothFR.traj,trajStruct1(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct1(j).allZSmoothFR,binWidth);               
        numTrajRemaining = numTraj - numSample;
        sampInd2 = randsample(numTrajRemaining,numSample);
        remainingInd = setdiff(1:numTraj,sampInd1);
        sampInd2 = remainingInd(sampInd2);
        trajStruct2(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd2);
        [trajStruct2(j).avgZSmoothFR.traj,trajStruct2(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct2(j).allZSmoothFR,binWidth);   
    end
        
end