function [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct,varargin)

% Get number/list of postures, targets, channels, and conditions
    dataType = 'zSmoothFR';
    assignopts(who,varargin);
    
    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    if strcmpi(dataType,'zSmoothFR')
        numChannels = size(trajStruct(1).avgZSmoothFR.traj,2);
    elseif strcmpi(dataType,'smoothFR')
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    end
    numConditions = size(trajStruct,2);

end