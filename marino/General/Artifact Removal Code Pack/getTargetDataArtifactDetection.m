function [targetData] = getTargetDataArtifactDetection(trialName,TrialTargets,trialRawData)

%% Get Trial Data
names = TrialTargets.names;
window = TrialTargets.window;

%% Match to case
foundCaseFlag = 0;
switch trialName
    %Delayed Center Out With Rewards
    case 'Delayed Center Out with Active Punishment'
        foundCaseFlag = 1;
        centerName = 'reachstart';
        targetName = 'reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    case 'Center Out 20160505'
        foundCaseFlag = 1;
        centerName = 'reachstart';
        targetName = 'reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    case 'Delayed Center Out Active Punishment Training'
        foundCaseFlag = 1;
        centerName = 'reachstart';
        targetName = 'reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
    %AutoMonkey
    case 'CenterOut_TouchBar_AutoMonkey'
        foundCaseFlag = 1;
        centerName = 'start';
        targetName = 't_bc_reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));

    %Center Out BC
    case 'CenterOut_TouchBar_BC'
        foundCaseFlag = 1;
        centerName = 'start';
        targetName = 't_bc_reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    case 'CenterOut_ForceBar_BC'
        foundCaseFlag = 1;
        centerName = 'start';
        targetName = 't_bc_reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    case 'CentetOut_BC_TouchBar'
        foundCaseFlag = 1;
        centerName = 'start';
        targetName = 't_start_grid_reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
    %Center Out Center BC - This case is unique, because I ran trial files
    %with the same trial name, but different target names.  So, check to
    %see whether the target names such as 't1_grid' are there.  If so,
    %those are the targets to use. Otherwise, use 't_start_grid_reach'
    case 'CenterOutCenter_BC_TouchBar'
        foundCaseFlag = 1;
        %Check to see whether "t1_grid," etc. are included
        possibleEndTargetNames = {'t1_grid','t2_grid','t3_grid','t4_grid','t5_grid','t6_grid','t7_grid','t8_grid'};
        containsOneOfTheseTargets = 0;
        for i=1:size(possibleEndTargetNames,2)
            possibleEndTargetName = possibleEndTargetNames{1,i};
            containsOneOfTheseTargets = containsOneOfTheseTargets + sum(cellfun(@(x) strcmpi(x,possibleEndTargetName),names));
            if containsOneOfTheseTargets > 0
                break
            end
        end
        %If so, that target is the end target.  Otherwise, it's
        %'t_start_grid_reach'
        if containsOneOfTheseTargets > 0
            targetName = possibleEndTargetName;
        else
            targetName = 't_start_grid_reach';
        end
        centerName = 'start';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
    %HC Center Out
    case 'HC_CenterOut'
        foundCaseFlag = 1;
        centerName = 't_center_touchbar';
        targetName = 't_reach_touchbar';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    
    %HC Grid Task
    case 'HC_GridTask'
        foundCaseFlag = 2;
        target1Name = 't_start_tube_small';
        target2Name = 't_end_tube_small';
        centerLocName = 'CoordinateOffset';
        AllPossibleTargets = trialRawData.Definitions.AllPossibleTargets;
        target1Ind  = min(find([cellfun(@(x) strcmpi(x,target1Name),names)]==1));
        target2Ind = min(find([cellfun(@(x) strcmpi(x,target2Name),names)]==1));
        centerLocInd = min(find([cellfun(@(x) strcmpi(x,centerLocName),AllPossibleTargets.names)]==1));
        %Assign Target Locations
        targetData.target1Loc = window(target1Ind,1:3);
        targetData.target1Size = window(target1Ind,4);
        targetData.target2Loc = window(target2Ind,1:3);
        targetData.target2Size = window(target2Ind,4);
        targetData.centerLoc = [AllPossibleTargets.xValues(centerLocInd,1),AllPossibleTargets.yValues(centerLocInd,1),AllPossibleTargets.zValues(centerLocInd,1)];
        % Flip Y
        targetData.target1Loc(1,2) = -targetData.target1Loc(1,2);
        targetData.target2Loc(1,2) = -targetData.target2Loc(1,2);
        targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
        %Get Target ID
        target2Vec = (targetData.target2Loc - targetData.centerLoc)';
        target2Vec(3,1) = 0;
        target2Ang = atan2d(norm(cross([1;0;0],target2Vec)),dot([1;0;0],target2Vec));
        if target2Vec(2,1) < 0
            target2Ang = 360-target2Ang;
        end
        target2Num = round(target2Ang./45) + 1;
        targetData.target2ID = target2Num; 
        
    %HC Grid Task w Tubes
    case 'HC_GridTask_Tube'
        foundCaseFlag = 2;
        target1Name = 't_start_tube_small';
        target2Name = 't_end_tube_small';
        centerLocName = 'CoordinateOffset';
        AllPossibleTargets = trialRawData.Definitions.AllPossibleTargets;
        target1Ind  = min(find([cellfun(@(x) strcmpi(x,target1Name),names)]==1));
        target2Ind = min(find([cellfun(@(x) strcmpi(x,target2Name),names)]==1));
        centerLocInd = min(find([cellfun(@(x) strcmpi(x,centerLocName),AllPossibleTargets.names)]==1));
        %Assign Target Locations
        targetData.target1Loc = window(target1Ind,1:3);
        targetData.target1Size = window(target1Ind,4);
        targetData.target2Loc = window(target2Ind,1:3);
        targetData.target2Size = window(target2Ind,4);
        targetData.centerLoc = [AllPossibleTargets.xValues(centerLocInd,1),AllPossibleTargets.yValues(centerLocInd,1),AllPossibleTargets.zValues(centerLocInd,1)];
        % Flip Y
        targetData.target1Loc(1,2) = -targetData.target1Loc(1,2);
        targetData.target2Loc(1,2) = -targetData.target2Loc(1,2);
        targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
                %Get Target ID
        target2Vec = (targetData.target2Loc - targetData.centerLoc)';
        target2Vec(3,1) = 0;
        target2Ang = atan2d(norm(cross([1;0;0],target2Vec)),dot([1;0;0],target2Vec));
        if target2Vec(2,1) < 0
            target2Ang = 360-target2Ang;
        end
        target2Num = round(target2Ang./45) + 1;
        targetData.target2ID = target2Num; 
        
    %Isometric Force 1D
    case 'IsometricForce_1D'
        foundCaseFlag = 1;
        centerName = 'forcecenter';
        targetName = 'forcetarget';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
    %GridTask Force Bar
    case 'GridTask_BC_ForceBar'
        foundCaseFlag = 2;
        centerName = 'start';
        target1Name = 't_start_grid_reach';
        target2Name = 't_end_grid_reach';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        target1Ind = min(find([cellfun(@(x) strcmpi(x,target1Name),names)]==1));
        target2Ind = min(find([cellfun(@(x) strcmpi(x,target2Name),names)]==1));
        % Assign Target Locations
        targetData.centerLoc = window(centerInd,1:3);
        targetData.centerSize = window(centerInd,4);
        targetData.target1Loc = window(target1Ind,1:3);
        targetData.target1Size = window(target1Ind,4);
        targetData.target2Loc = window(target2Ind,1:3);
        targetData.target2Size = window(target2Ind,4);
        % Flip Y
        targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
        targetData.target1Loc(1,2) = -targetData.target1Loc(1,2);
        targetData.target2Loc(1,2) = -targetData.target2Loc(1,2);
        %Get Target ID
        target1Vec = (targetData.target1Loc - targetData.centerLoc)';
        target1Vec(3,1) = 0;
        target1Ang = atan2d(norm(cross([1;0;0],target1Vec)),dot([1;0;0],target1Vec));
        if target1Vec(2,1) < 0
            target1Ang = 360-target1Ang;
        end
        target1Num = round(target1Ang./45) + 1;
        targetData.target1ID = target1Num; 
        target2Vec = (targetData.target2Loc - targetData.centerLoc)';
        target2Vec(3,1) = 0;
        target2Ang = atan2d(norm(cross([1;0;0],target2Vec)),dot([1;0;0],target2Vec));
        if target2Vec(2,1) < 0
            target2Ang = 360-target2Ang;
        end
        target2Num = round(target2Ang./45) + 1;
        targetData.target2ID = target2Num; 
end


%% Get TargetID
if foundCaseFlag == 1
    % Assign Target Locations
    targetData.centerLoc = window(centerInd,1:3);
    targetData.centerSize = window(centerInd,4);
    targetData.targetLoc = window(targetInd,1:3);
    targetData.targetSize = window(targetInd,4);
    targetData.workspaceCenter = targetData.centerLoc;
    
    %Center Target Locations
    targetData.targetLoc = targetData.targetLoc-targetData.workspaceCenter;
    targetData.centerLoc = targetData.centerLoc-targetData.workspaceCenter;
    
    % Flip Y
    targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
    targetData.targetLoc(1,2) = -targetData.targetLoc(1,2);
    targetData.workspaceCenter(1,2) = -1.*targetData.workspaceCenter(1,2);
    
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