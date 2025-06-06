function [trajStruct1,trajStruct2] = splitDataforCVKinematicComparison(trajStruct,minNumCondTrials,binWidth)       

    %Split trajStruct into two groups for cross-validation
    numSample = floor(minNumCondTrials/2); %Number of trajectories in each draw
    trajStruct1 = trajStruct;
    trajStruct2 = trajStruct;
    for j = 1:size(trajStruct,2)
        numTraj = size(trajStruct(j).allZSmoothFR,2);
        sampInd1 = randsample(numTraj,numSample);     
        trajStruct1(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd1);
        [trajStruct1(j).avgZSmoothFR.traj,trajStruct1(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct1(j),'zSmoothFR',binWidth);
        if isfield(trajStruct, 'allBciCursorTraj')
            trajStruct1(j).allBciCursorTraj = trajStruct(j).allBciCursorTraj(sampInd1);
            [trajStruct1(j).avgBciCursorTraj.traj,trajStruct1(j).avgBciCursorTraj.timestamps,trajStruct1(j).avgBciCursorTraj.CI95] = getAvgTraj(trajStruct1(j),'BciCursorTraj',binWidth); 
        end
        if isfield(trajStruct, 'allMarkerPos')
            trajStruct1(j).allMarkerPos = trajStruct(j).allMarkerPos(sampInd1);
            [trajStruct1(j).avgMarkerPos.traj,trajStruct1(j).avgMarkerPos.timestamps,trajStruct1(j).avgMarkerPos.CI95] = getAvgTraj(trajStruct1(j),'MarkerPos',binWidth); 
        end
        
        if isfield(trajStruct, 'allCursorPos')
            trajStruct1(j).allCursorPos = trajStruct(j).allCursorPos(sampInd1);
            [trajStruct1(j).avgCursorPos.traj,trajStruct1(j).avgCursorPos.timestamps,trajStruct1(j).avgCursorPos.CI95] = getAvgTraj(trajStruct1(j),'CursorPos',binWidth); 
        end
        
        if isfield(trajStruct, 'allForce')
            trajStruct1(j).allForce = trajStruct(j).allForce(sampInd1);
            [trajStruct1(j).avgForce.traj,trajStruct1(j).avgForce.timestamps,trajStruct1(j).avgForce.CI95] = getAvgTraj(trajStruct1(j),'Force',binWidth); 
        end
        
        
        numTrajRemaining = numTraj - numSample;
        sampInd2 = randsample(numTrajRemaining,numSample);
        remainingInd = setdiff(1:numTraj,sampInd1);
        sampInd2 = remainingInd(sampInd2);
        trajStruct2(j).allZSmoothFR = trajStruct(j).allZSmoothFR(sampInd2);
        [trajStruct2(j).avgZSmoothFR.traj,trajStruct2(j).avgZSmoothFR.timestamps] = getAvgTraj(trajStruct2(j),'zSmoothFR',binWidth);  
        if isfield(trajStruct, 'allBciCursorTraj')
            trajStruct2(j).allBciCursorTraj = trajStruct(j).allBciCursorTraj(sampInd2);
            [trajStruct2(j).avgBciCursorTraj.traj,trajStruct2(j).avgBciCursorTraj.timestamps,trajStruct2(j).avgBciCursorTraj.CI95] = getAvgTraj(trajStruct2(j),'BciCursorTraj',binWidth); 
        end
        if isfield(trajStruct, 'allMarkerPos')
            trajStruct2(j).allMarkerPos = trajStruct(j).allMarkerPos(sampInd2);
            [trajStruct2(j).avgMarkerPos.traj,trajStruct2(j).avgMarkerPos.timestamps,trajStruct2(j).avgMarkerPos.CI95] = getAvgTraj(trajStruct2(j),'MarkerPos',binWidth); 
        end
        if isfield(trajStruct, 'allCursorPos')
            trajStruct2(j).allCursorPos = trajStruct(j).allCursorPos(sampInd2);
            [trajStruct2(j).avgCursorPos.traj,trajStruct2(j).avgCursorPos.timestamps,trajStruct2(j).avgCursorPos.CI95] = getAvgTraj(trajStruct2(j),'CursorPos',binWidth); 
        end
        
        if isfield(trajStruct, 'allForce')
            trajStruct2(j).allForce = trajStruct(j).allForce(sampInd2);
            [trajStruct2(j).avgForce.traj,trajStruct2(j).avgForce.timestamps,trajStruct2(j).avgForce.CI95] = getAvgTraj(trajStruct2(j),'Force',binWidth); 
        end
                
    end   
end