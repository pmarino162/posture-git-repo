function [statesTraj,statesTimestamps] = getStatesTraj20211210(trialData,trialInclStates,dataType,binWidth,varargin)

%Variable Arguments 
timeRelToTrialStart = false;
chMeans = [];
chStd = [];
rmCh = [];
assignopts(who,varargin);

%Determine trialType
trialName = trialData.trialName;
numTrialTypes = size(trialInclStates,2);
for trialTypeInd = 1:numTrialTypes
   if strcmpi(trialInclStates(trialTypeInd).trialName,trialName)
       trialType = trialTypeInd;
   end
end

%Create empty outputs and variables
statesTraj = [];
statesTimestamps = [];
startTime = [];
endTime = [];

%Get State Info
stateTransitions = double(trialData.stateData.stateTransitions);
stateNames = trialData.stateData.stateNames;

%Get kinData (if it's there)
if isfield(trialData,'kinData')
   kinData = trialData.kinData; 
end

%Get start and end times
inclStates = trialInclStates(trialType).inclStates;
for event =1:2
    eventType = inclStates{event}{1};
    eventName = inclStates{event}{2};
    eventOccurence = inclStates{event}{3};
    eventAddTime = inclStates{event}{4};
    
    %Get event time
    if strcmpi(eventType,'state')
        if strcmpi(eventName,'Last')
            eventTime = double(stateTransitions(2,end));
        else
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
        end
    elseif strcmpi(eventType,'kin')
        eventTime = getfield(kinData,eventName);
    end
    
    %Add time to event time; store align time (time of event 1)
    if ~isempty(eventTime)
        if event == 1
            alignTime = eventTime;
            startTime = eventTime + eventAddTime;
            if startTime > alignTime
                alignTime = startTime;
            end
            if startTime < 0 
               error('trial startTime < 0') 
            end
        elseif event == 2 
            endTime = eventTime + eventAddTime;
        end
    end
    
end

%Get trialStartTime and trialEndTime
trialStartTime = double(stateTransitions(2,1));
trialEndTime = double(stateTransitions(2,end));

%Extract state traj and state timesteps if you got a valid start and end
%time
if ~isempty(startTime) && ~isempty(endTime)
    switch dataType
        
        %For smoothed neural data, bin @ 1ms, then smooth, then downsample
        %(if necessary)
        case {'smoothFR','zSmoothFR'}
        tempBinWidth = 1;
        %Design Gaussian kernel for convolution (25ms std dev)
        sigma = 50/tempBinWidth; %sigma in number of bins 
        windowSideLength = ceil(3*sigma); %half width of gaussian kernal in bins
        kernel = normpdf([-windowSideLength:1:windowSideLength],0,sigma);
        kernel = kernel.*(1/sum(kernel));
        %Adjust start and end times if there isn't enough data to smooth
        if startTime < trialStartTime + (windowSideLength+0.5)*tempBinWidth
            startTime = trialStartTime + (windowSideLength+0.5)*tempBinWidth;
            if startTime > alignTime
                alignTime = startTime;
            end
        end
        if endTime > trialEndTime - (windowSideLength+0.5)*tempBinWidth
            endTime = trialEndTime - (windowSideLength+0.5)*tempBinWidth;
        end
        %Bin @ 1ms (including enough extra bins to properly smooth)
        numBackBins = ceil((abs(alignTime-startTime)/tempBinWidth)-0.5) + windowSideLength;
        numForwardBins = ceil(((endTime-alignTime)/tempBinWidth)-0.5) + windowSideLength;
        timestamps = -1*(numBackBins*tempBinWidth):tempBinWidth:numForwardBins*tempBinWidth;
        timestamps = timestamps + alignTime;
        binEdges = [timestamps-tempBinWidth/2,timestamps(end)+tempBinWidth/2];
        binnedSpikes = getBinnedSpikes20211210(trialData.spikes,binEdges);
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
        if binWidth ~= 1
            if mod(binWidth,1) ~= 0
               error('Error: binWidth for smoothFR must be integer for downsampling') 
            end
            alignInd = find(timestamps==alignTime);
            keepInd = unique([alignInd:-binWidth:1,alignInd:binWidth:length(timestamps)]);
            statesTraj = binnedSpikes(keepInd,:);
            statesTimestamps = timestamps(keepInd);
        else
            statesTraj = binnedSpikes;
            statesTimestamps = timestamps;
        end
        
        %Z-score smoothedFR
        if strcmpi(dataType,'zSmoothFR')  
            if isempty(chMeans) || isempty(chStd)
                error('Must specify z-score parameters')
            end
            statesTraj = (statesTraj-chMeans)./chStd;
            statesTraj(:,rmCh) = [];
        end
        
        %For binned neural data, just bin activity
        case 'binFR'
        numBackBins = ceil((abs(alignTime-startTime)/binWidth)-0.5);
        numForwardBins = ceil(((endTime-alignTime)/binWidth)-0.5);
        timestamps = -1*(numBackBins*binWidth):binWidth:numForwardBins*binWidth;
        timestamps = timestamps + alignTime;
        binEdges = [timestamps-binWidth/2,timestamps(end)+binWidth/2];
        statesTraj = getBinnedSpikes20211210(trialData.spikes,binEdges);
        statesTraj = statesTraj./(binWidth/1000); %convert to Hz
        statesTimestamps = timestamps;
        
        %For non-neural data streams, snip, then downsample
        case{'marker','markerVel','force','forceCursor','eyePos','pupil'}
        if strcmpi(dataType,'force')
           timestamps = trialData.force.time;
           traj = trialData.force.force;
        elseif strcmpi(dataType,'forceCursor')
           timestamps = trialData.force.time;
           traj = trialData.force.forceCursor;
        elseif strcmpi(dataType,'marker')
           timestamps = trialData.marker.time;
           traj = trialData.marker.position;
        elseif strcmpi(dataType,'eyePos')
           timestamps = trialData.eye.time;
           traj = trialData.eye.position;
        elseif strcmpi(dataType,'pupil')
           timestamps = trialData.eye.time;
           traj = trialData.eye.pupil;
        end
        timestampsMask = (timestamps >= startTime & timestamps <= endTime);
        statesTraj = traj(timestampsMask,:);
        statesTimestamps = double(timestamps(timestampsMask));
        %Downsample (if binWidth ~= 1)
        if binWidth ~= 1
            if mod(binWidth,1) ~= 0
               error('Error: binWidth must be integer for downsampling') 
            end
            alignInd = find(statesTimestamps==alignTime);
            keepInd = unique([alignInd:-binWidth:1,alignInd:binWidth:length(statesTimestamps)]);
            statesTraj = statesTraj(keepInd,:);
            statesTimestamps = statesTimestamps(keepInd);
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

end