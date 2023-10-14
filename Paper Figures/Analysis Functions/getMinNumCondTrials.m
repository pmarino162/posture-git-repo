function [minNumCondTrials] = getMinNumCondTrials(trajStruct)    
    
        %Get minimum number of trials in any condition
        numCondTraj = nan(1,size(trajStruct,2));
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allZSmoothFR,2);
           numCondTraj(i) = numTraj;
        end
        [minNumCondTrials,~] = min(numCondTraj);
        
end