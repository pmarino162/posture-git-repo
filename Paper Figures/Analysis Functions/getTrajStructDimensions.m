function [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct)

% Get number/list of postures, targets, channels, and conditions

    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgZSmoothFR.traj,2);
    numConditions = size(trajStruct,2);

end