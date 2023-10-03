function [statesTraj,statesTimestamps] = getStatesTraj(trialData,trialInclStates,dataType,binWidth,kernelStdDev,varargin)

% This function takes in data from one trial and returns a trajectory  
% of a user-specified type during a user-specified window of the trial 

%% Variable Arguments 
    timeRelToTrialStart = false;
    zScoreParams = [];
    rmSort = [];
    assignopts(who,varargin);

%% Determine trialType
    trialName = trialData.trialName;
    numTrialTypes = size(trialInclStates,2);
    for trialTypeInd = 1:numTrialTypes
       if strcmpi(trialInclStates(trialTypeInd).trialName,trialName)
           trialType = trialTypeInd;
       end
    end

%% Create empty outputs
    statesTraj = []; statesTimestamps = [];

%% Get State Info
    stateTransitions = double(trialData.stateData.stateTransitions);
    stateNames = trialData.stateData.stateNames;

%% Get kinData (if it's there)
    if isfield(trialData,'kinData')
       kinData = trialData.kinData; 
    end

%% Get start and end times
    inclStates = trialInclStates(trialType).inclStates;
    trialTime = trialData.time;
    for event = 1:2
        eventType = inclStates{event}{1};
        eventName = inclStates{event}{2};
        eventOccurence = inclStates{event}{3};
        eventAddTime = inclStates{event}{4};

        %Get event time
        if strcmpi(eventType,'state')
            stateInd = max(find([cellfun(@(x) strcmpi(x,eventName),stateNames)]==1));          
            if ismember(stateInd,stateTransitions(1,:))
                if strcmpi(eventOccurence,'first')
                    eventTime = double(stateTransitions(2,min(find(stateTransitions(1,:)==stateInd))));
                elseif strcmpi(eventOccurence,'last')
                    eventTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd))));
                end
            else
                error('Did not find event name in state table')
            end
        elseif strcmpi(eventType,'kin')
            eventTime = getfield(kinData,eventName);
        end

        %Add time to event time; store align time (time of event 1)
        if ~isempty(eventTime)
            if event == 1
                alignTime = eventTime;
                startTime = eventTime + eventAddTime;
                if startTime < trialTime(1)
                   error('startTime < first trial timestamp') 
                end
            elseif event == 2 
                endTime = eventTime + eventAddTime;
                if endTime > trialTime(end)
                   error('endTime > last trial timestamp') 
                end
            end
        end

    end

%% Get trialStartTime and trialEndTime
    trialStartTime = trialTime(1);
    trialEndTime = trialTime(end);

%% Extract state traj and state timesteps if you got a valid start and end time
    if ~isempty(startTime) && ~isempty(endTime)
        switch dataType               
            %Smooths neural data (binned at 1ms by default) then downsamples to desired 'binWidth'
            case {'smoothFR','zSmoothFR'}
                %Design Gaussian kernel for convolution
                tempBinWidth = 1; %For readability of calculations below. This just reflects that neural data is binned at 1ms before processing
                sigma = kernelStdDev/tempBinWidth; %sigma in number of bins (kernelStdDev and tempBinWidth are in ms)
                windowSideLength = ceil(3*sigma); %half width of gaussian kernal in bins
                kernel = normpdf([-windowSideLength:1:windowSideLength],0,sigma);
                kernel = kernel.*(1/sum(kernel));
                %Adjust start and end times if there isn't enough data to smooth
                if startTime < trialStartTime + windowSideLength*tempBinWidth
                    startTime = trialStartTime + windowSideLength*tempBinWidth;
                end
                if endTime > trialEndTime - windowSideLength*tempBinWidth
                    endTime = trialEndTime - windowSideLength*tempBinWidth;
                end
                %Get appropriate time window of binnedSpikes
                timestamps = startTime-windowSideLength:endTime+windowSideLength;
                trialTimeMask = trialTime>=timestamps(1) & trialTime<=timestamps(end);
                binnedSpikes = trialData.spikes(trialTimeMask,:);
                %Convolve with Guassian kernel
                for sort = 1:size(binnedSpikes,2)
                    binnedSort = binnedSpikes(:,sort);
                    temp = nan(size(binnedSort,1),size(binnedSort,2));
                    temp(windowSideLength+1:length(temp)-windowSideLength,1) = conv(binnedSort,kernel,'valid');
                    binnedSpikes(:,sort) = temp;
                end
                %Remove extra bins and timestamps
                rmInd = timestamps<startTime | timestamps>endTime;
                timestamps(rmInd) = [];
                binnedSpikes(rmInd,:) = [];
                binnedSpikes = binnedSpikes./(tempBinWidth/1000); %Convert to Hz
                %Downsample (if binWidth ~= 1)
                statesTraj = binnedSpikes;
                statesTimestamps = timestamps;
                if binWidth ~= 1
                    [statesTraj,statesTimestamps] = downsample(statesTraj,statesTimestamps,binWidth,alignTime);
                end
                %Z-score smoothedFR
                if strcmpi(dataType,'zSmoothFR')  
                    %Error handling
                    if isempty(zScoreParams)
                        error('Must specify z-score parameters')
                    end
                    tempZScoreParams = zScoreParams([zScoreParams.binWidth]==binWidth & [zScoreParams.kernelStdDev]==kernelStdDev);
                    if isempty(tempZScoreParams)
                        error('Did not pass in z-score parameters for current binWidth and kernelStdDev')
                    end
                    %Implementation
                    sortMean = tempZScoreParams.sortMean;
                    sortStd = tempZScoreParams.sortStd;
                    statesTraj = (statesTraj-sortMean)./sortStd;
                    statesTraj(:,rmSort) = [];
                end

            %Didn't end up using unsmoothed binned FRs for this study
            case 'binFR'

            %For non-neural data streams, snip to appropriate time window, then downsample
            case{'marker','markerPos','markerVel','force','forceCursor','bciCursorTraj','bciCursorPos','bciCursorVel'}
                if strcmpi(dataType,'force')
                   timestamps = trialData.force.time;
                   traj = trialData.force.force;
                elseif strcmpi(dataType,'forceCursor')
                   timestamps = trialData.force.time;
                   traj = trialData.force.forceCursor;
                elseif strcmpi(dataType,'marker')
                   timestamps = trialData.marker.time;
                   traj = trialData.marker.position;
                elseif strcmpi(dataType,'markerPos')
                   timestamps = trialData.marker.time;
                   traj = trialData.marker.position;
                elseif strcmpi(dataType,'markerVel')
                   timestamps = trialData.marker.time;
                   traj = trialData.marker.velocity;
                elseif strcmpi(dataType,'bciCursorTraj')
                   timestamps = trialData.Decoder.timestamps;
                   traj = trialData.Decoder.cursorTraj;
                elseif strcmpi(dataType,'bciCursorPos')
                   timestamps = trialData.Decoder.posTime;
                   traj = trialData.Decoder.position;
                elseif strcmpi(dataType,'bciCursorVel')
                   timestamps = trialData.Decoder.velTime;
                   traj = trialData.Decoder.velocity;
                end
                timestampsMask = (timestamps >= startTime & timestamps <= endTime);
                statesTraj = traj(timestampsMask,:);
                statesTimestamps = double(timestamps(timestampsMask));
                if ~isempty(statesTraj)
                    %Downsample (if binWidth ~= 1) and save
                    if binWidth ~= 1
                        [statesTraj,statesTimestamps] = downsample(statesTraj,statesTimestamps,binWidth,alignTime);
                    end 
                end
        end
    end

    %If ~timeRelToStart, adjust timestamps so that they're relative to align
    %time
    if ~isempty(statesTimestamps)
        if ~timeRelToTrialStart
            statesTimestamps = statesTimestamps-alignTime;
        end
    end

%% Local function for downsampling statesTraj and statesTimestamps
 function [statesTraj,statesTimestamps] = downsample(statesTraj,statesTimestamps,binWidth,alignTime)
    if mod(binWidth,1) ~= 0
       error('Error: binWidth must be integer for downsampling') 
    end
    keepTimestamps = unique([statesTimestamps(1):binWidth:alignTime,alignTime,alignTime:binWidth:statesTimestamps(end)]);
    keepInd = ismember(statesTimestamps,keepTimestamps);
    statesTraj = statesTraj(keepInd,:);
    statesTimestamps = statesTimestamps(keepInd);
 end  
 
end