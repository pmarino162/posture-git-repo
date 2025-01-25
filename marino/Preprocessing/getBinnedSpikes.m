function [binTimes,channelSpikes,allChannelSpikeBins,allSortSpikeBins] = getBinnedSpikes(binWidth,spikes,stateTransitions,channelSpikes,getSorts)
   
    %This function bins spikes for each sort on each channel.  It also
    %combines bins for all sorts on each channel.  It uses the histcounts
    %function for binning, so some extra steps have to be
    %taken to work around this function's limitations
        
    %% Get channel and sort info 
    channelList = [channelSpikes.channel];
    uniqueChannelList = unique(channelList);
    sortList = nan(1,size(channelSpikes,2));
    for i = 1:size(channelSpikes,2)
        if isnumeric(channelSpikes(i).sort)
            sortList(i) = channelSpikes(i).sort;
        end
    end
    numChannels = size(uniqueChannelList,2); 
    numValidSorts = sum(~isnan(sortList) & sortList~= 31 & sortList ~= 0);
    
    %% Set up bins
    startTime = double(stateTransitions(2,1));
    endTime = double(stateTransitions(2,end));
    binEdges = startTime:binWidth:endTime;
    binTimes = uint32(binEdges(2:end));
    numBins = size(binTimes,2);
    allChannelSpikeBins = zeros(numBins,numChannels,'uint16');
    allSortSpikeBins = zeros(numBins,numValidSorts,'uint16');
    %Before running, add extra last bin so that histcounts won't include both
    %edges of final bin
    binEdges = [binEdges,binEdges(end)];
    %If you want to include the right bin edge, add 1 to each bin edge
    %(this works since timestamps are always integers).  Otherwise, bins
    %will be left-edged
    binEdges = binEdges + 1;
    
    %% For each channel, bin spikes for each sort and for all sorts 
    if ~isempty(spikes)
        channelSpikeBinsInd = 1;
        sortSpikeBinsInd = 1;
        for channel = uniqueChannelList
            %Get channel data
            chSpikes = spikes([spikes.channel]==channel);
            totalChSpikes = zeros(1,numBins,'uint16');
            %Bin spikes for each sort.  Count up channel total
            for sortInd = 1:size(chSpikes,2)
               sortID = chSpikes(sortInd).sort;
               timestamps = double(chSpikes(sortInd).timestamps);
               binnedSpikeMat = uint16(histcounts(timestamps,binEdges));
               %Delete extra last bin
               binnedSpikeMat(end) = [];
               %Add sort data to channelSpikes
               if getSorts == true
                   channelSpikesInd = find(channelList==channel & sortList==sortID);
                   channelSpikes(channelSpikesInd).binnedSpikes = binnedSpikeMat;
                   %Save to allSortSpikeBins
                   if sortID == 0 || sortID == 31
                       %If sort 0 or sort 31, don't save to allSortSpikeBins
                   else
                       allSortSpikeBins(:,sortSpikeBinsInd) = binnedSpikeMat';
                       sortSpikeBinsInd = sortSpikeBinsInd + 1;
                   end
               end               
               totalChSpikes = totalChSpikes + binnedSpikeMat;
            end
            %Save channel totals
            if getSorts == true
                channelSpikes(channelSpikesInd+1).binnedSpikes = totalChSpikes;
            else
                channelSpikesInd = find(channelList==channel);
                channelSpikes(channelSpikesInd).binnedSpikes = totalChSpikes;
            end
            allChannelSpikeBins(:,channelSpikeBinsInd) = totalChSpikes';
            channelSpikeBinsInd = channelSpikeBinsInd + 1;
            %Save total channel spikes to allChannelSpikeBins and sorted spikes to
            %allSortSpikeBins if the channel was turned on (if it was off,
            %these matrices will still have column(s) of zeros for the
            %channel)
% 
%                 %Spike Counts
%                 
%                 for sortInd = 1:size(chSpikes,2)
%                    allSortSpikeBins(:,sortSpikeBinsInd) = channelSpikes(channel).sorts(sortInd).bins';
%                    sortSpikeBinsInd = sortSpikeBinsInd + 1;
%                 end
%             end
%         end
         
        end
        
    else
    %If there are no spikes, don't fill channelSpikes or allSpikeBins
    end
end