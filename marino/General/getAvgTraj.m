function [avgTraj,timestamps] = getAvgTraj(trajStruct,extTrajStruct,method)

%% Determine timestamps for avgTraj; preallocate
    numCondTrials = size(trajStruct,2);
    
    %Get start times, end times, and sample period
    startTimes = nan(1,numCondTrials);
    endTimes = nan(1,numCondTrials);
    allSamplePeriods = [];
    for trial = 1:numCondTrials
        trialTimestamps = trajStruct(trial).timestamps;
        if ~isempty(trialTimestamps)
            startTimes(trial) = trialTimestamps(1);
            endTimes(trial) = trialTimestamps(end);
        end
        allSamplePeriods = [allSamplePeriods,diff(trialTimestamps)];
    end
    
    %Get median start and end times
    startTime = median(startTimes,'omitnan');
    endTime = median(endTimes,'omitnan');
    
    %Get mode inter-timestamp length 
    samplePeriod = mode(allSamplePeriods);
    
    %Create timestamp vector for average trajectory 
        numForwardBins = ceil((endTime-startTime-samplePeriod/2)/samplePeriod);
        timestamps = startTime:samplePeriod:startTime + numForwardBins*samplePeriod;
        %old method for centering at 0
% %     numBackBins = ceil((abs(startTime)/samplePeriod)-1/2);
% %     numForwardBins = ceil((endTime/samplePeriod)-0.5);
% %     timestamps = -1*(numBackBins*samplePeriod):samplePeriod:numForwardBins*samplePeriod;
        binEdges = [timestamps-samplePeriod/2,timestamps(end)+samplePeriod/2];
    
    %Preallocate average trajectory
    firstValidTrajInd = find(endTimes>0,1,'first');
    traj = trajStruct(firstValidTrajInd).traj;
    numDims = size(traj,2);
    avgTraj = zeros(size(timestamps,2),numDims);
    
%% Take Average
    for i = 1:size(binEdges,2)-1
        binStart = binEdges(i);
        binEnd = binEdges(i+1);
        numPts = 0;
        for trial = 1:numCondTrials    
            trialTimestamps = extTrajStruct(trial).timestamps;
            traj = extTrajStruct(trial).traj;
            %Interpolation Method
            if strcmpi(method,'interp')
%                 avgTraj(i,:) = interp1(binStart,traj,) + avgTraj(i,:);
%                 numPts = numPts + 1;
            %Binning Method
            elseif strcmpi(method,'binning')
                trajInd = find(trialTimestamps >= binStart & trialTimestamps < binEnd);
                %If there are multiple samples from a trial within a bin, take
                %average 
                if size(trajInd,2) > 1
                    avgTraj(i,:) = mean(traj(trajInd,:)) + avgTraj(i,:);
                    numPts = numPts + 1;
                elseif size(trajInd,2) == 1
                    avgTraj(i,:) = traj(trajInd,:) + avgTraj(i,:);
                    numPts = numPts + 1;
                end
            end
        end
        avgTraj(i,:) = avgTraj(i,:)/numPts;
    end


end