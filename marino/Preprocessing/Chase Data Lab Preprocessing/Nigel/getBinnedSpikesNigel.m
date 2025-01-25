function [binTimes,channelSpikes,allChannelSpikeBins,allSortSpikeBins] = getBinnedSpikesNigel(binWidth,spikes,stateTransitions,channelList,channelSpikes,edge)
   
    %This function bins spikes for each sort on each channel.  It also
    %combines bins for all sorts on each channel.  
    
    %Chase Lab Data includes spike times which are non-integers. In order
    %to retain the ability to specify bin edgedness (i.e., right or left)
    %with this type of data, I am not using histcounts
    
    %% Nigel Array Info
    numPossibleChannels = 96;
    numPossibleSortsPerChannel = 4;
    
    %% Set up bins
    startTime = double(stateTransitions(2,1));
    endTime = double(stateTransitions(2,end));
    binEdges = startTime:binWidth:(endTime+binWidth);
    if strcmpi(edge,'left')
        binTimes = uint32(binEdges(1:end-1));
    elseif  strcmpi(edge,'right')
        binTimes = uint32(binEdges(2:end));
    end
    numBins = size(binTimes,2);
    
    %% For each channel that was used, bin spikes for each sort and for all sorts 
    if ~isempty(spikes)
        % Preallocate allChannelSpikeBins and allSortSpikeBins
        allChannelSpikeBins = zeros(numBins,numPossibleChannels,'uint16');
        allSortSpikeBins = zeros(numBins,numPossibleChannels*numPossibleSortsPerChannel,'uint16');
        channelSpikesInd = 1;
        for channel = channelList
            %Get channel data
            chSpikes = spikes([spikes.channel]==channel);
            totalChSpikes = zeros(1,numBins,'uint16');
            numSorts = size(chSpikes,2);
            %Bin spikes for each sort.  Count up channel total
            for sortInd = 1:numSorts
               sortID = chSpikes(sortInd).sort;
               timestamps = double(chSpikes(sortInd).timestamps);
               binnedSpikeMat = zeros(1,numBins);
               %Left-edged bins
               if strcmpi(edge,'left')
                   for i = 1:size(timestamps,1)
                      bin = floor(timestamps(i)/binWidth)+1; 
                      binnedSpikeMat(bin) = binnedSpikeMat(bin) + 1;
                   end
               %Right-edged bins
               elseif  strcmpi(edge,'right')
                   for i = 1:size(timestamps,1)
                      bin = ceil(timestamps(i)/binWidth); 
                      binnedSpikeMat(bin) = binnedSpikeMat(bin) + 1;
                   end
               end
               %Convert to 16 bit
               binnedSpikeMat = uint16(binnedSpikeMat);
               %Fill Struct
               channelSpikes(channelSpikesInd).sorts(sortInd).sortID = uint16(chSpikes(sortInd).sort);
               channelSpikes(channelSpikesInd).sorts(sortInd).bins = binnedSpikeMat;
               totalChSpikes = totalChSpikes + binnedSpikeMat;
            end
            channelSpikes(channelSpikesInd).allSorts = totalChSpikes;
            %Spike Counts
            allChannelSpikeBins(:,channel) = totalChSpikes';
            for sortInd = 1:numSorts
               allSortSpikeBins(:,4*(channel-1)+sortInd) = channelSpikes(channelSpikesInd).sorts(sortInd).bins';
            end
            channelSpikesInd = channelSpikesInd + 1;
        end
    else
         %If there are no spikes, don't fill channelSpikes or allChannelSpikeBins
         allChannelSpikeBins = [];
         allSortSpikeBins = [];
    end
end