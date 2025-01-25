function [binnedSpikes] = getBinnedSpikes20211210(spikes,binEdges)

%Preallocate binnedSpikes
numBins = length(binEdges) - 1;
numSorts = size(spikes,2);
binnedSpikes = zeros(numBins,numSorts);

%Bin spikes (right-edged)
for sort = 1:numSorts
    timestamps = spikes(sort).timestamps;
    for bin = 1:numBins
       binnedSpikes(bin,sort) = sum(timestamps>binEdges(bin) & timestamps<=binEdges(bin+1));
    end
end

end