%This function takes in a data matrix from a trial and returns the portion of that
%matrix which was recorded during the specified occurrences of the specified
%inclStates

function [statesTraj,statesTimestamps,startTime] = getStatesTraj(trialData,trialInclStates,dataType,varargin)
    %Variable Arguments 
    timeRelToTrialStart = false;
    assignopts(who,varargin);
    
    %Get relevant data for dataType 
    switch dataType
        case 'marker'
            X = double(trialData.marker.position);
            timestamps = double(trialData.marker.time);
        case 'markerVel'
            X = double(trialData.marker.velocity);
            timestamps = double(trialData.marker.time);
        case 'force'
            X = double(trialData.force.force);
            timestamps = double(trialData.force.time);
        case 'absoluteForce'
            X = double(trialData.force.absoluteForcesAndTorques);
            timestamps = double(trialData.force.time);
        case 'forceCursor'
            X = double(trialData.force.forceCursor);
            timestamps = double(trialData.force.time);   
        case 'eye'
            X1 = double(trialData.eye.position);
            X2 = double(trialData.eye.pupil);
            timestamps = double(trialData.eye.time); 
        case {'allChannelSpikeBins','allChannelSmoothedFR'}
            X = double(trialData.spikes.allChannelSpikeBins);
            timestamps = double(trialData.spikes.binTimes);
        case {'allSortSpikeBins','allSortSmoothedFR'}
            X = double(trialData.spikes.allSortSpikeBins);
            timestamps = double(trialData.spikes.binTimes);
        case 'GPFA'
            X = trialData.spikes.GPFA;
            timestamps = double(trialData.spikes.GPFABinTimes);
        case 'GPFAVel'
            X = trialData.spikes.GPFAVel;
            timestamps = double(trialData.spikes.GPFAVelTimes);
        case 'GPFAAllSpikeBins'
            X = trialData.spikes.GPFAAllSpikeBins;
            timestamps = double(trialData.spikes.GPFABinTimes);    
        case 'PCProj'
            X = trialData.spikes.PCProj;
            timestamps = double(trialData.spikes.binTimes);
        case 'PCVel'
            X = trialData.spikes.PCVel;
            timestamps = double(trialData.spikes.binTimes(1:end-1));
        case 'Decoder'
            X1 = trialData.Decoder.cursorTraj;
            X2 = trialData.Decoder.GPFATraj;
            X3 = trialData.Decoder.WorkTraj;
            X4 = trialData.Decoder.NullTraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'decoderGPFATraj'
            X = trialData.Decoder.GPFATraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'decoderCursorTraj'
            X = trialData.Decoder.cursorTraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'decoderWorkTraj'
            X = trialData.Decoder.WorkTraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'N00WorkTraj'
            X = trialData.Decoder.N00WorkTraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'I30WorkTraj'
            X = trialData.Decoder.I30WorkTraj;
            timestamps = double(trialData.Decoder.timestamps);
        case 'E30WorkTraj'
            X = trialData.Decoder.E30WorkTraj;
            timestamps = double(trialData.Decoder.timestamps);      
        case 'simDecoder'
            X = trialData.simDecoder.simCursorTraj;
            timestamps = double(trialData.simDecoder.timestamps);
        case 'decoderSpikeBins'
            X = double(trialData.Decoder.rawSpikeBins);
            timestamps = double(trialData.Decoder.timestamps);
    end
    
    %Determine trialType
    trialName = trialData.trialName;
    numTrialTypes = size(trialInclStates,2);
    for trialTypeInd = 1:numTrialTypes
       if strcmpi(trialInclStates(trialTypeInd).trialName,trialName)
           trialType = trialTypeInd;
       end
    end
        
    %Get State Info
    inclStates = trialInclStates(trialType).inclStates;
    numInclStates = size(inclStates,2);
    stateTransitions = double(trialData.stateData.stateTransitions);
    stateNames = trialData.stateData.stateNames;

    %Get kinData
    
    
    %For each inclState, find start and end time (if state occurred). Update timestampsMask
    %accordingly
    timestampsMask = false(1,length(timestamps));
    for state = 1:numInclStates
        stateInd = max(find([cellfun(@(x) strcmpi(x,inclStates{1,state}),stateNames)]==1));
        occurrence = trialInclStates(trialType).inclOccurrence{1,state};
        if ismember(stateInd,stateTransitions(1,:))
            %Get State Start and End Times 
                if strcmpi(occurrence,'first')
                    stateStartTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==stateInd))));
                    stateEndTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==stateInd))+1));
                    startTime = stateStartTime;
                    endTime = stateEndTime;
                elseif strcmpi(occurrence,'last')
                    stateStartTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd))));
                    stateEndTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd))+1));
                    startTime = stateStartTime;
                    endTime = stateEndTime;
                end
            %Manually add time to these if specified by user
                if isfield(trialInclStates,'addTimeToBeginning')
                    addTimeToBeginning = trialInclStates(trialType).addTimeToBeginning{1,state};
                    if size(addTimeToBeginning,2) == 1
                        startTime = stateStartTime + addTimeToBeginning(1,1);
                    elseif size(addTimeToBeginning,2) == 2
                        if isfield(trialInclStates,'addTimeToEnd') && trialInclStates(trialType).addTimeToEnd{1,state}~=0
                            error('If addTimeToBeginning is 2x2 for a state, addTimeToEnd must equal 0 for that state')
                        end
                        endTime = stateStartTime + addTimeToBeginning(1,2);
                        startTime = stateStartTime + addTimeToBeginning(1,1); 
                    else
                        error('addTimeToBeginning must contain a 1x1 or 1x2 vector for each state')
                    end
                        
                end
                if isfield(trialInclStates,'addTimeToEnd')
                    addTimeToEnd = trialInclStates(trialType).addTimeToEnd{1,state};
                    if size(addTimeToEnd,2) == 1
                        endTime = stateEndTime + addTimeToEnd(1,1);
                    elseif size(addTimeToEnd,2) == 2
                        if isfield(trialInclStates,'addTimeToBeginning') && trialInclStates(trialType).addTimeToBeginning{1,state}~=0
                            error('If addTimeToEnd is 2x2 for a state, addTimeToBeginning must equal 0 for that state')
                        end
                        startTime = stateEndTime + addTimeToEnd(1,1);
                        endTime = stateEndTime + addTimeToEnd(1,2);
                    else
                        error('addTimeToEnd must contain a 1x1 or 1x2 vector for each state')                      
                    end
                end
            %Update timestampsMask    
                timestampsMask = (timestamps >= startTime & timestamps <= endTime) | timestampsMask;
            %Save start time for first state
                if state == 1
                   firstStateStartTime = stateStartTime; 
                end
        end
    end

    %Get Trajectory
    switch dataType
        %For most dataTypes, snip out the relevant range
        case {'marker','markerVel','force','absoluteForce','forceCursor','allChannelSpikeBins','allSortSpikeBins','GPFA','GPFAVel',...
                'GPFAAllSpikeBins','PCProj','PCVel','simDecoder','decoderSpikeBins','decoderCursorTraj',...
                'decoderWorkTraj','decoderGPFATraj','N00WorkTraj','I30WorkTraj','E30WorkTraj'}
            statesTraj = X(timestampsMask,:);
        case 'eye'
            statesTraj.position = X1(timestampsMask,:);
            statesTraj.pupil = X2(timestampsMask,:);
        case 'Decoder'
            statesTraj.cursorTraj = X1(timestampsMask,:);
            statesTraj.GPFATraj = X2(timestampsMask,:);
            statesTraj.WorkTraj = X3(timestampsMask,:);
            statesTraj.NullTraj = X4(timestampsMask,:);
        %For smoothed firing rates, smooth spike bins, then snip
        case {'allChannelSmoothedFR','allSortSmoothedFR'}
            binWidth = mode(diff(timestamps));
            X = X./(binWidth/1000); %convert to Hz
            %Convolve w Gaussian (25ms std dev)
            sigma = 25/binWidth; %sigma in number of bins 
            windowSideLength = ceil(3*sigma); %half width of gaussian kernal in bins
            kernel = normpdf([-windowSideLength:1:windowSideLength],0,sigma);
            kernel = kernel.*(1/sum(kernel));
            for sort = 1:size(X,2)
                Xsort = X(:,sort);
                temp = nan(size(Xsort,1),size(Xsort,2));
                temp(windowSideLength+1:length(temp)-windowSideLength,1) = conv(Xsort,kernel,'valid');
                X(:,sort) = temp;
            end
            statesTraj = X(timestampsMask,:);
    end
    statesTimestamps = double(timestamps(timestampsMask));
    if ~isempty(statesTimestamps)
        if ~timeRelToTrialStart
            statesTimestamps = statesTimestamps-firstStateStartTime;
        end
    end
end