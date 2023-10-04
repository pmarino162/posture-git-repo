function [minNumCondTrials] = getMinNumCondTrials(trajStruct)    

        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allZSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end
        [minNumCondTrials,i] = min(numCondTraj);
        
end