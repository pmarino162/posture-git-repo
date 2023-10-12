function [totalNumTrials] = getTotalNumTrials(trajStruct,varargin)    

    totalNumTrials = 0;
    for i = 1:size(trajStruct,2)
        if isfield(trajStruct,'allZSmoothFR')
            numTraj = size(trajStruct(i).allZSmoothFR,2);
        elseif isfield(trajStruct,'allSmoothFR')
            numTraj = size(trajStruct(i).allSmoothFR,2);
        end
       totalNumTrials = totalNumTrials + numTraj;
    end
        
end