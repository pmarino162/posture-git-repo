function [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct)

% Get number/list of postures, targets, channels, and conditions
    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    if isfield(trajStruct,'avgZSmoothFR')
        numChannels = size(trajStruct(1).avgZSmoothFR.traj,2);
    elseif isfield(trajStruct,'avgSmoothFR')
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    elseif isfield(trajStruct,'avgSingleBinFR')
        numChannels = size(trajStruct(1).avgSingleBinFR.traj,2);
    end
    numConditions = size(trajStruct,2);

end