%This function preallocates the structure I use to hold binned spiking data
% exclZero: if true, exclude channels that only have 1 sort whose ID is 0
% getSorts: if true, bin spikes for each sort separately. if false, combine
% sorts for each channel
function [spikes] = preallocateSpikesStruct(spikes,exclZero,getSorts)


        channelSpikes = struct('channel',uint16([]),'sort',uint16([]),'binnedSpikes',uint16([]));
        
        channelList = unique([spikes.channel]);
        %Exclude channels
        if exclZero == true
            zeroCh = [];
            for channel = channelList
                chSpikes = spikes([spikes.channel]==channel);
                chSorts = [chSpikes.sort];
                if size(chSorts,2) == 1 && chSorts(1) == 0
                    zeroCh = [zeroCh,channel];
                end
            end
            channelList = setdiff(channelList,zeroCh);
        end
        channelList = setdiff(channelList,exclCh);
        structInd = 1;
        %Get Sorts 
        if getSorts == true
            for channel = channelList
                chSpikes = spikes([spikes.channel]==channel);
                chSorts = [chSpikes.sort];
                for sort = chSorts
                    channelSpikes(structInd).channel = channel;
                    channelSpikes(structInd).sort = sort;
                    structInd = structInd + 1;
                end
                channelSpikes(structInd).channel = channel;
                channelSpikes(structInd).sort = 'all';
                structInd = structInd + 1;
            end
            sortList = rmfield(channelSpikes,'binnedSpikes');
            sortCell = {sortList.sort};
            rowsToDelete = cellfun(@(x) strcmpi(x,'all'),sortCell);
            sortList(rowsToDelete) = [];
            sortList([sortList.sort]==0 | [sortList.sort]==31) = [];
            spikes = struct('binTimes',uint32([]),'channelSpikes',channelSpikes,'allSortSpikeBins',uint16([]),'sortList',sortList);
        %Get Only Channels
        else
            for channel = channelList
                channelSpikes(structInd).channel = channel;
                channelSpikes(structInd).sort = 'all';
                structInd = structInd + 1;
            end
            sortList = rmfield(channelSpikes,'binnedSpikes');
            spikes = struct('binTimes',uint32([]),'channelSpikes',channelSpikes,'allChannelSpikeBins',uint16([]),'sortList',sortList);
        end



end