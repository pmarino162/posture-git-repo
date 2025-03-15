function [maxNumCondTrials] = getMaxNumCondTrials(trajStruct)    
    
        %Get minimum number of trials in any condition
        numCondTraj = nan(1,size(trajStruct,2));
        if isfield(trajStruct,'allZSmoothFR')
            for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allZSmoothFR,2);
               numCondTraj(i) = numTraj;
            end
        
        elseif isfield(trajStruct,'allSmoothFR')
            for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               numCondTraj(i) = numTraj;
            end
            
        end
        [maxNumCondTrials,~] = max(numCondTraj);
        
end