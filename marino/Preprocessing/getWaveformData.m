function [waveforms] = getWaveformData(spikes)

%This function organizes the waveform data from the spikes struct into a
%format that I prefer.


if ~isempty(spikes)
    %Get unique list of channels
    channelList = unique([spikes.channel]);
    for channel = channelList
        %Get channel data
        chSpikes = spikes([spikes.channel]==channel);
        numSorts = size(chSpikes,2);
        %For each sort, collect waveforms
        allTimestamps = [];
        allWaveforms = [];
        for sortInd = 1:numSorts
           %Get Data
           timestamps = chSpikes(sortInd).timestamps;
           allTimestamps = [allTimestamps,timestamps];
           spikeWaveforms = chSpikes(sortInd).spikeWaveforms;
           allWaveforms = [allWaveforms,spikeWaveforms];
           %Fill struct
           waveforms(channel).sorts(sortInd).sortID = chSpikes(sortInd).sort;
           waveforms(channel).sorts(sortInd).timestamps = timestamps;
           waveforms(channel).sorts(sortInd).waveforms = spikeWaveforms;
        end
        %Sort the waveforms from all sorts by time
        [waveforms(channel).allSorts.timestamps,I] = sort(allTimestamps);
        waveforms(channel).allSorts.waveforms = allWaveforms(:,I);
    end
else
    %If there are no spikes, return an empty array
    waveforms = [];
end


end