function [totalNumTrials] = getTotalNumTrials(trajStruct)    

        totalNumTrials = 0;
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allZSmoothFR,2);
           totalNumTrials = totalNumTrials + numTraj;
        end
        
end