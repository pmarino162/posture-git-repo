function [minNumTimestamps] = getMinNumTimestamps(trajStruct)    

        % Get minimum number of timestamps in condition averages
        numTimestamps = nan(1,size(trajStruct,2));
        for i = 1:size(trajStruct,2)          
            if isfield(trajStruct,'avgZSmoothFR')
                numTimestamps(i) = length(trajStruct(i).avgZSmoothFR.timestamps);
            elseif isfield(trajStruct,'avgSmoothFR')
                numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
            end
        end
        [minNumTimestamps,~] = min(numTimestamps);
        
end