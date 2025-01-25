function [Data] = standardizeSpikesStruct(Data)
%This function gets all possible unit/sort combinations from an experiment,
%then updates the spikes struct on every trial to reflect all of those
%units/sorts

%% Get all possible sorts
    spikes = Data(1).spikes;
    allPossibleSorts = [[spikes.channel]',[spikes.sort]'];  
    numTrials = size(Data,2);
    for trial = 2:numTrials
        spikes = Data(trial).spikes;
        trialSorts = [[spikes.channel]',[spikes.sort]'];  
        newSorts = setdiff(trialSorts,allPossibleSorts,'rows');
        allPossibleSorts = vertcat(allPossibleSorts,newSorts);
        allPossibleSorts = sortrows(allPossibleSorts);
    end
    
%% Create new empty spikes struct
    spikes = struct('channel',[],'sort',[],'timestamps',zeros(1,1000));
    numSorts = size(allPossibleSorts,1);
    emptySorts = true(1,numSorts);
    spikes = repmat(spikes,1,numSorts);
    for i = 1:numSorts
       spikes(i).channel = allPossibleSorts(i,1);
       spikes(i).sort = allPossibleSorts(i,2);
    end
    
%% Fill it on every trial
%     for trial = 1:size(Data,2)
%         trial
%         trialSpikes = Data(trial).spikes;
%         tempSpikes = spikes;
%         for i = 1:numSorts
%            channel = tempSpikes(i).channel;
%            sort = tempSpikes(i).sort;
%            %If trial contains data for that channel and sort, add the
%            %timestamps to tempSpikes
%            if any([trialSpikes.channel]==channel & [trialSpikes.sort]==sort)
%                 timestamps = trialSpikes([trialSpikes.channel]==channel & [trialSpikes.sort]==sort).timestamps;
%                 tempSpikes(i).timestamps = timestamps;
%            else
%                tempSpikes(i).timestamps = [];
%            end
%         end
%         %Replace trial's spikes struct with tempSpikes
%         Data(trial).spikes = tempSpikes;
%     end

    for trial = 1:numTrials
        trial
        trialSpikes = Data(trial).spikes;
        tempSpikes = spikes;
        
        %Get indices of trialSpikes that have non-empty timestmaps
        hasTimestamps = zeros(1,size(trialSpikes,2));
        for i = 1:size(trialSpikes,2)
            if ~isempty(trialSpikes(i).timestamps)
               hasTimestamps(i) = 1; 
            end
        end 
        hasTimestampsInd = find(hasTimestamps);
        
        %For each row of trialSpikes with timestamps, transfer to
        %tempSpikes. Convert all others to empty.
        tempSpikesRowsChanged = zeros(1,size(hasTimestampsInd,2));
        tempSpikesRowChangedInd = 1;
        for i = hasTimestampsInd
            channel = trialSpikes(i).channel;
            sort = trialSpikes(i).sort;
            timestamps = trialSpikes(i).timestamps;
            tempSpikesRow = find([tempSpikes.channel]==channel & [tempSpikes.sort]==sort);
            tempSpikesRowsChanged(1,tempSpikesRowChangedInd) = tempSpikesRow;
            tempSpikes(tempSpikesRow).timestamps = timestamps;
            tempSpikesRowChangedInd = tempSpikesRowChangedInd + 1;
        end
        %Keep track of which sorts aren't empty
        emptySorts(tempSpikesRowsChanged) = false;
        %Clear out placeholders for sorts with no timestamps on this trial
        for j = setdiff(1:numel(tempSpikes),tempSpikesRowsChanged)
            tempSpikes(j).timestamps = [];
        end
        %Replace trial's spikes struct with tempSpikes
        Data(trial).spikes = tempSpikes;
    end
    
    %Clear Empty Sorts
    for trial = 1:numTrials
        Data(trial).spikes = Data(trial).spikes(~emptySorts);
    end
    
end