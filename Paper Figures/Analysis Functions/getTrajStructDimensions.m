function [postureList,numPostures,targetList,numTargets,numChannels,numConditions,varargout] = getTrajStructDimensions(trajStruct)

% Get number/list of postures, targets, channels, and conditions
    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    if isfield(trajStruct,'task')
       taskList = unique([trajStruct.task]);
       numTasks = size(taskList,2);
       varargout(1) = taskList;
       varargout(2) = numTasks;
    end
    if isfield(trajStruct,'avgZSmoothFR')
        numChannels = size(trajStruct(1).avgZSmoothFR.traj,2);
    elseif isfield(trajStruct,'avgSmoothFR')
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    elseif isfield(trajStruct,'avgSingleBinFR')
        numChannels = size(trajStruct(1).avgSingleBinFR.traj,2);
    end
    numConditions = size(trajStruct,2);

end