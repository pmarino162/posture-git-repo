function [procData,P] = preprocessBatistaData(Data, varargin)

%Preprocesses raw data from Batista lab. Uses my old function,
%preprocessDataMarino to do most of this, then does some extra stuff. 

%% Assign optional arguments 
    checkPhasespaceSync = false;
    droppedMarkerTrials = 'all';
    checkDecode = false;
    removeScreenFreezeTrials = false;
    removeBadPhasespace = true;
    assignopts(who,varargin);

%% Do most of the work with preprocessDataMarino
    procData = preprocessDataMarino(Data,'CHECKPHASESPACESYNC',checkPhasespaceSync,'DROPPEDMARKERTRIALS',droppedMarkerTrials,'CHECKDECODE',checkDecode,'REMOVESCREENFREEZETRIALS',removeScreenFreezeTrials);

%% Remove duplicate phasespace samples and trials for which phasespace wasn't saved properly
    %if there are duplicate psp samples, remove them
    numTrials = size(procData,2);
    for trial = 1:numTrials
        rowsToDelete = [];
        rawPositions = procData(trial).TrialData.Marker.rawPositions;
        if ~isempty(rawPositions)
            time = rawPositions(:,6);
            numRows = size(rawPositions,1);
            for row = 1:numRows-1
                if time(row+1,1) == time(row,1)
                    rowsToDelete = [rowsToDelete,row+1];
                end
            end
            rawPositions(rowsToDelete,:) = [];
            procData(trial).TrialData.Marker.rawPositions = rawPositions;
        end
    end

    %remove any trials for which the phasespace data wasn't saved properly
    if removeBadPhasespace
        numTrials = size(procData,2);
        rmTrials = [];
        for trial = 1:numTrials
            rawPositions = procData(trial).TrialData.Marker.rawPositions;
            if ~isempty(rawPositions)
                pspEndTime = rawPositions(end,6);
                trialEndTime = procData(trial).TrialData.stateTransitions(2,end);
                if trialEndTime > pspEndTime + 500
                    rmTrials = [rmTrials,trial];
                end
            end
        end
        procData(rmTrials) = [];
    end
    
end