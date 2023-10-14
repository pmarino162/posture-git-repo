function [trajStruct1,trajStruct2] = splitDataforCV(trajStruct)       

    %Split trajStruct into two groups for cross-validation
    numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
    trajStruct1 = trajStruct;
    trajStruct2 = trajStruct;
    for j = 1:size(trajStruct,2)
        numTraj = size(trajStruct(j).allSmoothFR,2);
        sampInd1 = randsample(numTraj,numSample);
        trajStruct1(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd1);
        [trajStruct1(j).avgSmoothFR.traj,trajStruct1(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct1(j).allSmoothFR,binWidth);               
        numTrajRemaining = numTraj - numSample;
        sampInd2 = randsample(numTrajRemaining,numSample);
        remainingInd = setdiff(1:numTraj,sampInd1);
        sampInd2 = remainingInd(sampInd2);
        trajStruct2(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd2);
        [trajStruct2(j).avgSmoothFR.traj,trajStruct2(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct2(j).allSmoothFR,binWidth);   
    end
        
end