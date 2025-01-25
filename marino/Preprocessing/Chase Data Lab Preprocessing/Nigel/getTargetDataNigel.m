function [targetData] = getTargetDataNigel(trialName,header,trialRawData)

%foundCaseFlag == 0 indicates that a case was not found for the trial type.
%foundCaseFlag == 1 indicates that a case was found for the trial type, and 
%   leads to a set of instructions that is common to many trial types. 
%foundCaseFlag == 2 indicates that a case was found for the trial type, but
%   that further instructions will be specific to that trial type.


%% Match to case
foundCaseFlag = 0;
switch trialName
    %Delayed Center Out With Rewards
    case 'Nigel Posture BC Center Out'
        foundCaseFlag = 1;
        % Assign Target Locations
        targetData.centerLoc = trialRawData.StartPos;
        targetData.centerSize = header.cout.TargetRadius;
        targetData.targetLoc = trialRawData.TargetPos;
        targetData.targetSize = header.cout.TargetRadius;
        targetData.workspaceCenter = targetData.centerLoc;
        %Center Target Locations
        targetData.targetLoc = targetData.targetLoc-targetData.workspaceCenter;
        targetData.centerLoc = targetData.centerLoc-targetData.workspaceCenter;
end


%% Get TargetID
if foundCaseFlag == 1
    %Get Target ID
    targetVec = (targetData.targetLoc - targetData.centerLoc)';
    targetVec(3,1) = 0;
    targetAng = atan2d(norm(cross([1;0;0],targetVec)),dot([1;0;0],targetVec));
    if targetVec(2,1) < 0
        targetAng = 360-targetAng;
    end
    targetNum = round(targetAng./45) + 1;
    targetData.targetID = targetNum; 
    
elseif foundCaseFlag == 0
    fprintf('Didn''t find trial case for targets\n')
    targetData = [];
end

end