function [totalNumTrials] = getTotalNumTrials(trajStruct,varargin)    
    dataType = 'zSmoothFR';
    assignopts(who,varargin);

        totalNumTrials = 0;
        for i = 1:size(trajStruct,2)
           
           
            if strcmpi(dataType,'zSmoothFR')
                numTraj = size(trajStruct(i).allZSmoothFR,2);
            elseif strcmpi(dataType,'smoothFR')
                numTraj = size(trajStruct(i).allSmoothFR,2);
            end
           totalNumTrials = totalNumTrials + numTraj;
        end
        
end