function [targetData] = getTargetData(trialName,TrialTargets,trialRawData)

%foundCaseFlag == 0 indicates that a case was not found for the trial type.
%foundCaseFlag == 1 indicates that a case was found for the trial type, and 
%   leads to a set of instructions that is common to many trial types. 
%foundCaseFlag == 2 indicates that a case was found for the trial type, but
%   that further instructions will be specific to that trial type.

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
    case 'CenterOut_20181112' %Prez Choking
        foundCaseFlag = 1;
        centerName = 'start';
        targetName = 'reachtarget';
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
    case {'CenterOut_TouchBar_BC','BCI Center Out'}
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
    %with the same trial name but different target names.  So, check to
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
        dateStr = num2str(trialRawData.Overview.date);   
        if strcmpi(dateStr(1:8),'20190729') %On this date, center target name shows up as 'No Target Name'
            centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1 | [cellfun(@(x) strcmpi(x,'neural decoder'),names)]==1));
        else
            centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        end       
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
    %HC Center Out
    case 'HC_CenterOut'
        foundCaseFlag = 1;
        centerName = 't_center_touchbar';
        targetName = 't_reach_touchbar';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    
    %Earl Multi-Posture Grid Reaching, Delayed Center Out, and Posture Device Movements 
    case {'GridReaching','GridReachingCatch','Delayed Center Out 20210621','Delayed Center Out Catch 20210621','Posture Device Active Movements',...
            'DelayedCenterOut20210828'}
        foundCaseFlag = 1;
        centerName = 'center';
        targetName = 'target';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
    
    case 'HC_CenterOut_ForceBar_20200314'
        foundCaseFlag = 2;
        targetName = 'hc_targets';
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        %Get Coordinate Offset for Current Posture
        AllPossibleTargets = trialRawData.Definitions.AllPossibleTargets;
        AllPossibleTargetNames = AllPossibleTargets.names;
        coordinateOffsetInd = min(find([cellfun(@(x) strcmpi(x,'CoordinateOffset'),AllPossibleTargetNames)]==1));
        coordinateOffset = [AllPossibleTargets.xValues(coordinateOffsetInd,1), AllPossibleTargets.yValues(coordinateOffsetInd,1), AllPossibleTargets.zValues(coordinateOffsetInd,1)];
        % Assign Target Locations
        targetData.targetLoc = window(targetInd,1:3);
        targetData.targetSize = window(targetInd,4);
        targetData.centerLoc = coordinateOffset;
        targetData.workspaceCenter = [150 350 -680]; %Makes all phasespace data relative to neutral workspace center
        targetData.coordinateOffset = coordinateOffset;
        %Center Target Locations
        targetData.targetLoc = targetData.targetLoc-targetData.workspaceCenter;
        targetData.centerLoc = targetData.centerLoc-targetData.workspaceCenter + [0 -15 0];
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
        
    %Isometric Force 1D - In Jan 2020, I ran some MVC experiments using
    %trials with the same trial name. These experiments used a target named
    %'forcetargetmvc.'  If this target was used, make it the targetName.
    case 'IsometricForce_1D'
        foundCaseFlag = 2;
        %Check to see whether 'forcetargetmvc' was included
        possibleEndTargetNames = {'forcetargetmvc'};
        containsOneOfTheseTargets = 0;
        for i=1:size(possibleEndTargetNames,2)
            possibleEndTargetName = possibleEndTargetNames{1,i};
            containsOneOfTheseTargets = containsOneOfTheseTargets + sum(cellfun(@(x) strcmpi(x,possibleEndTargetName),names));
            if containsOneOfTheseTargets > 0
                break
            end
        end
        %If so, that target is the end target.  Otherwise, it's
        %'forcetarget'
        if containsOneOfTheseTargets > 0
            targetName = possibleEndTargetName;
        else
            targetName = 'forcetarget';
        end
        centerName = 'forcecenter';
        centerInd  = min(find([cellfun(@(x) strcmpi(x,centerName),names)]==1));
        targetInd = min(find([cellfun(@(x) strcmpi(x,targetName),names)]==1));
        
        % Assign Target Locations
        targetData.centerLoc = window(centerInd,1:3);
        targetData.centerSize = window(centerInd,4:5);
        targetData.targetLoc = window(targetInd,1:3);
        targetData.targetSize = window(targetInd,4:5);
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
        
    %GridTask Force Bar
    case {'GridTask_BC_ForceBar','GridTask_BC_ForceBar_16_Target','GridTask_CO_Across_BC_ForceBar'}
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
        targetData.workspaceCenter = targetData.centerLoc;
        targetData.target1Loc = window(target1Ind,1:3);
        targetData.target1Size = window(target1Ind,4);
        targetData.target2Loc = window(target2Ind,1:3);
        targetData.target2Size = window(target2Ind,4);
        % Flip Y
        targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
        targetData.target1Loc(1,2) = -targetData.target1Loc(1,2);
        targetData.target2Loc(1,2) = -targetData.target2Loc(1,2);
        targetData.workspaceCenter(1,2) = -1.*targetData.workspaceCenter(1,2);
        %Get Target ID
        target1Vec = (targetData.target1Loc - targetData.centerLoc)';
        target1Vec(3,1) = 0;
        target1Ang = atan2d(norm(cross([1;0;0],target1Vec)),dot([1;0;0],target1Vec));
        if target1Vec(2,1) < 0
            target1Ang = 360-target1Ang;
        end
        target2Vec = (targetData.target2Loc - targetData.centerLoc)';
        target2Vec(3,1) = 0;
        target2Ang = atan2d(norm(cross([1;0;0],target2Vec)),dot([1;0;0],target2Vec));
        if target2Vec(2,1) < 0
            target2Ang = 360-target2Ang;
        end
        if strcmpi(trialName,'GridTask_BC_ForceBar') || strcmpi(trialName,'GridTask_CO_Across_BC_ForceBar')
            target1Num = round(target1Ang./45) + 1;
            target2Num = round(target2Ang./45) + 1;
        elseif strcmpi(trialName,'GridTask_BC_ForceBar_16_Target')
            target1Num = round(target1Ang./22.5) + 1;
            target2Num = round(target2Ang./22.5) + 1;
        end
        targetData.target1ID = target1Num; 
        targetData.target2ID = target2Num; 
        
    %HC Grid Task
    case 'HC_MemoryGuided_TargetRing_ForceBar'
        foundCaseFlag = 2;
        bcTarget1Name = 'bc_ta';
        bcTarget2Name = 'bc_tb';
        hcTarget1Name = 'hc_t1';
%         centerLocName = 'CoordinateOffset';
        AllPossibleTargets = trialRawData.Definitions.AllPossibleTargets;
        bcTarget1Ind  = min(find([cellfun(@(x) strcmpi(x,bcTarget1Name),names)]==1));
        bcTarget2Ind  = min(find([cellfun(@(x) strcmpi(x,bcTarget2Name),names)]==1));
        hcTarget1Ind  = min(find([cellfun(@(x) strcmpi(x,hcTarget1Name),names)]==1));
%         centerLocInd = min(find([cellfun(@(x) strcmpi(x,centerLocName),AllPossibleTargets.names)]==1));
        %Assign Target Locations
        targetData.bcTarget1Loc = window(bcTarget1Ind,1:3);
        targetData.bcTarget1Size = window(bcTarget1Ind,4);
        targetData.bcTarget2Loc = window(bcTarget2Ind,1:3);
        targetData.bcTarget2Size = window(bcTarget2Ind,4);
        targetData.hcTarget1Loc = window(hcTarget1Ind,1:3);
        targetData.hcTarget1Size = window(hcTarget1Ind,4);
        targetData.centerLoc = [-55,330,0];
        targetData.workspaceCenter = targetData.centerLoc;        
        % Flip Y
        targetData.bcTarget1Loc(1,2) = -targetData.bcTarget1Loc(1,2);
        targetData.bcTarget2Loc(1,2) = -targetData.bcTarget2Loc(1,2);
        targetData.hcTarget1Loc(1,2) = -targetData.hcTarget1Loc(1,2);
        targetData.centerLoc(1,2) = -targetData.centerLoc(1,2);
        targetData.workspaceCenter(1,2) = -targetData.workspaceCenter(1,2);
        %Center Target Locations
        targetData.bcTarget1Loc = targetData.bcTarget1Loc-targetData.workspaceCenter;
        targetData.bcTarget2Loc = targetData.bcTarget2Loc-targetData.workspaceCenter;
        targetData.hcTarget1Loc = targetData.hcTarget1Loc-targetData.workspaceCenter;
        targetData.centerLoc = targetData.centerLoc-targetData.workspaceCenter;
    
        %Get Target ID
        hcTarget1Vec = (targetData.hcTarget1Loc - targetData.centerLoc)';
        hcTarget1Vec(3,1) = 0;
        hcTarget1Ang = atan2d(norm(cross([1;0;0],hcTarget1Vec)),dot([1;0;0],hcTarget1Vec));
        if hcTarget1Vec(2,1) < 0
            hcTarget1Ang = 360-hcTarget1Ang;
        end
        hcTarget1Num = round(hcTarget1Ang./45) + 1;
        targetData.hcTarget1ID = hcTarget1Num; 
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