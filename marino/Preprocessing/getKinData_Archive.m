function [kinData] = getKinData(trialName,trialRawData,trialData,centerMarker)


switch trialName
    %Earl Multi-Posture Grid Reaching, Multi-tasking reaching
    case {'GridReaching','GridReachingCatch','Delayed Center Out 20210621','Delayed Center Out Catch 20210621',...
            'DelayedCenterOut20210828','HC_CenterOut_ForceBar_20200314'}
        %Set up trial include states
        %trialData.trialNum
        trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {trialName};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
        %Get state, target, and interval data
        stateNames = trialRawData.Parameters.stateNames;
        stateTransitions = trialData.stateData.stateTransitions;
        intervals = [trialRawData.Parameters.StateTable.Interval];
        if strcmpi(trialName,'HC_CenterOut_ForceBar_20200314')
            %State
            centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Touch Bar Hold'),stateNames)]==1));
            targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Reach'),stateNames)]==1));
            targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Hold'),stateNames)]==1));
            %Interval
            centerHoldTime = intervals(centerHoldID).length;
            targetHoldTime = intervals(targetHoldID).length;
            %Target
            centerLoc = trialData.targetData.centerLoc;
        else
            %State
            centerHoldID = max(find([cellfun(@(x) strcmpi(x,'Center Hold'),stateNames)]==1));
            delayID = max(find([cellfun(@(x) strcmpi(x,'Delay'),stateNames)]==1));
            targetAcqID = max(find([cellfun(@(x) strcmpi(x,'Target Acquire'),stateNames)]==1));
            targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
            %Interval
            centerHoldTime = intervals(centerHoldID).length;
            delayLength = intervals(delayID).length;
            targetHoldTime = intervals(targetHoldID).length;
            %Target
            centerLoc = trialData.targetData.centerLoc;
            centerSize = trialData.targetData.centerSize;
        end

        %Get additional info for valid reaches (made it to hold)
        if any(ismember(targetHoldID,stateTransitions(1,:)))
            goTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetAcqID))));
            targetEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldID))));
            %Get phasespace data
            markerTime = trialData.marker.time;
            markerPos = trialData.marker.position;
            if ~centerMarker
                workspaceCenter = trialData.targetData.workspaceCenter;
                markerPos = markerPos - workspaceCenter(1,1:2);
            end
            markerDist = vecnorm(markerPos');
            %Get peak speed info from reach
            if strcmpi(trialName,'HC_CenterOut_ForceBar_20200314')
                trialInclStates(1).inclStates = {'Reach'};
            else
                trialInclStates(1).inclStates = {'Target Acquire'};
            end
            trialInclStates(1).inclOccurrence = {'first'};
            [markerVel,velTime,~] = getStatesTraj(trialData,trialInclStates,'markerVel','timeRelToTrialStart',true);
            speed = vecnorm(markerVel');
            [peakSpeed,peakSpeedInd] = max(speed);
            peakSpeedTime = velTime(peakSpeedInd);
            moveOnsetSpeed = 0.2*peakSpeed;
            %Walk backwards from peak speed to get movement onset index
            i = peakSpeedInd;
            while speed(i) > moveOnsetSpeed
                i = i-1;
            end
            moveOnsetInd = i;
            moveOnsetTime = velTime(moveOnsetInd);  
        
%                     figure
%                     plot(velTime,markerVel(:,1))
%                     hold on
%                     plot(velTime,markerVel(:,2))
%                     title('velocity')
%                     
%                     figure
%                     plot(velTime,speed)
%                     title('speed')
                    
%             %Get center exit time
%             slack = 5; % mm outside of center radius for reaction time 
%             if isempty(goTime)
%                 firstIndAfterCenterExit = [];
%                 lastIndBeforeCenterExit = [];
%             else
%                 firstIndAfterCenterExit = min(find(markerDist>centerSize+slack & markerTime>goTime));
%                 firstIndTime = markerTime(firstIndAfterCenterExit);
%                 if ~isempty(firstIndAfterCenterExit)
%                     lastIndBeforeCenterExit = max(find(markerDist<centerSize+slack & markerTime<firstIndTime));
%                 end
%             end
%             if isempty(firstIndAfterCenterExit)
%                 centerExitTime = [];
%             else
%                 centerExitTime = interp1(markerDist(lastIndBeforeCenterExit:firstIndAfterCenterExit),markerTime(lastIndBeforeCenterExit:firstIndAfterCenterExit),centerSize+slack);
%             end
        end
        
        %Store in kinData
        if strcmpi(trialName,'HC_CenterOut_ForceBar_20200314')
        else
            kinData.delayLength = delayLength;
            if any(ismember(targetHoldID,stateTransitions(1,:)))
                %kinData.centerExitTime = centerExitTime;
                %kinData.reactTimeCenterExit = centerExitTime-goTime;
            end
        end
        kinData.centerHoldTime = centerHoldTime;
        kinData.targetHoldTime = targetHoldTime;
        if any(ismember(targetHoldID,stateTransitions(1,:)))
            kinData.moveOnsetTime = moveOnsetTime;
            kinData.reactTimeMoveOnset = moveOnsetTime-goTime;
            kinData.acqTime =  targetEnterTime-goTime;
            kinData.reachTime = targetEnterTime-moveOnsetTime;
            kinData.peakSpeed = peakSpeed;
        end
        
    %Earl BCI Tasks
    case {'GridTask_CO_Across_BC_ForceBar','BCI Center Out'}
        %Set up trial include states
        trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
        trialInclStates(1).trialName = {trialName};
        trialInclStates(1).addTimeToBeginning = {0};
        trialInclStates(1).addTimeToEnd = {0};
        %Get state data
        stateNames = trialRawData.Parameters.stateNames;
        stateTransitions = trialData.stateData.stateTransitions;
        step1ID = max(find([cellfun(@(x) strcmpi(x,'Step 1'),stateNames)]==1));
        step2ID = max(find([cellfun(@(x) strcmpi(x,'Step 2'),stateNames)]==1));
        successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
        %Get state transitions and step acquisition times
        if ismember(step1ID,stateTransitions(1,:))
            step1Time = double(stateTransitions(2,min(find(stateTransitions(1,:)==step1ID))));
        end
        if ismember(step2ID,stateTransitions(1,:))
            step2Time = double(stateTransitions(2,min(find(stateTransitions(1,:)==step2ID))));
            step1AcqTime = step2Time - step1Time;
            kinData.step1AcqTime = step1AcqTime;
            kinData.movementTime = step1AcqTime;
        end
        if ismember(successID,stateTransitions(1,:)) 
            successTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==successID))));
            if ismember(step2ID,stateTransitions(1,:))
                step2AcqTime = successTime - step2Time;
                kinData.step2AcqTime = step2AcqTime;
            else  
                kinData.movementTime = successTime - step1Time;
            end
        end
        if ~ismember(step2ID,stateTransitions(1,:)) & ~ismember(successID,stateTransitions(1,:)) 
           kinData = []; 
        end
%         %Get movement onset, reaction, and movement times
%         if ismember(step2ID,stateTransitions(1,:))
%             trialInclStates(1).inclStates = {'Step 1'};
%             trialInclStates(1).inclOccurrence = {'first'};
%             centerLoc = trialData.targetData.centerLoc;
%             target1Loc = trialData.targetData.target1Loc;
%             target1Vec = target1Loc(1,1:2)-centerLoc(1,1:2);
%             target1Dist = norm(target1Vec);
%             [cursorTraj,timestamps,~] = getStatesTraj(trialData,trialInclStates,'decoderCursorTraj','timeRelToTrialStart',true);
%             target1Proj = dot(cursorTraj,repmat(target1Vec,length(timestamps),1),2)./target1Dist;
%             threshInd = min(find(target1Proj > 0.5*target1Dist));
%             if threshInd == 1
%                 moveOnsetTime = timestamps(1);
%             else
%                 moveOnsetTime = interp1(target1Proj(threshInd-1:threshInd),timestamps(threshInd-1:threshInd),0.5*target1Dist);
%             end
%             kinData.moveOnsetTime = round(moveOnsetTime);
%             kinData.rxnTime = moveOnsetTime - step1Time;
%             kinData.movementTime = step2Time-moveOnsetTime;
%         end

     
    case {'CenterOut_ForceBar_BC'}
        %Get state data
        stateNames = trialRawData.Parameters.stateNames;
        stateTransitions = trialData.stateData.stateTransitions;
        touchBarBCID = max(find([cellfun(@(x) strcmpi(x,'Touch Bar BC'),stateNames)]==1));
        successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
        %Get state transitions
        goTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==touchBarBCID))));
        targetEnterTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==successID))));
        acqTime = targetEnterTime-goTime;
        %Store in kinData
        kinData.acqTime = acqTime;
 
        
    %Isometric Force
    case{'IsometricForce_1D'}
    %Set up trial include states
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {trialName};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0};
    %Get state data
    stateNames = trialData.stateData.stateNames;
    stateTransitions = trialData.stateData.stateTransitions;
    targetID = max(find([cellfun(@(x) strcmpi(x,'Target'),stateNames)]==1));
    targetHoldID = max(find([cellfun(@(x) strcmpi(x,'Target Hold'),stateNames)]==1));
    successID = max(find([cellfun(@(x) strcmpi(x,'Success with Reward'),stateNames)]==1));
    %Get state transitions and step acquisition times
    if ismember(targetID,stateTransitions(1,:))
        targetTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetID))));
    end
    if ismember(targetHoldID,stateTransitions(1,:))
        targetHoldTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==targetHoldID))));
        targetAcqTime = targetHoldTime - targetTime;
        kinData.targetAcqTime = targetAcqTime;
    else
        kinData = [];
    end
    %Get movement onset, reaction, and movement times if trial made it to
    %target hold
    if ismember(targetHoldID,stateTransitions(1,:))
        %Get Y force trajectory after go cue
        trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
        [targetForceTraj,timestamps] = getStatesTraj20211210(trialData,trialInclStates,'force',1,'timeRelToTrialStart',true);
        YTargetForceTraj = targetForceTraj(:,2);
        %Walk backwards from target hold time to nearest time when force
        %went 20% towards end force above baseline mean in correct direction. Interpolate to improve estimate
        threshold = 0.2*(YTargetForceTraj(end)-YTargetForceTraj(1)) + YTargetForceTraj(1);
        dir = sign(YTargetForceTraj(end)-YTargetForceTraj(1));
        i = length(YTargetForceTraj);
        if dir > 0
            while YTargetForceTraj(i) > threshold
               i = i-1; 
            end
            moveOnsetTime = interp1(YTargetForceTraj(i:i+1),timestamps(i:i+1),threshold);
        elseif dir < 0 
            while YTargetForceTraj(i) < threshold
               i = i-1; 
            end
            moveOnsetTime = interp1(YTargetForceTraj(i:i+1),timestamps(i:i+1),threshold);
        end
        kinData.moveOnsetTime = round(moveOnsetTime);
        kinData.rxnTime = moveOnsetTime - targetTime;
    end
    
    otherwise
        kinData = [];
end

end