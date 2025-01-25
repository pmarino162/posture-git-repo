function [GPFATraj,WorkTraj,NullTraj] = getDecoderGPFATraj(Decoder)

    %Get Relevant Parameters
    W = Decoder.Parameters.W;
    NullSpace = null(W);
    M = Decoder.Parameters.M;
    CRinv = Decoder.Parameters.CRinv;
    c = Decoder.Parameters.c;
    d = Decoder.Parameters.d;
    
    %Demean and project spike bins
    rawSpikeBins = double(Decoder.rawSpikeBins);
    rawSpikeBinsT = rawSpikeBins';
    spikes_demean = rawSpikeBinsT-d;
    spikeProj = CRinv*spikes_demean;
    
    %Smooth projected spikes and save into GPFATraj
    GPFATraj = zeros(size(rawSpikeBinsT,2),10);
    smoothVec = zeros(70,1);
    for t = 1:size(rawSpikeBinsT,2)
        smoothVec(1:60) = smoothVec(11:70);
        smoothVec(61:70) = spikeProj(:,t);
        GPFATraj(t,:) = transpose(M*smoothVec);
    end
        
    %Use GPFATraj to get workspace and null space projection 
    WorkTraj = transpose(W*GPFATraj');
    numNullDims = size(NullSpace,2);
    NullTraj = zeros(size(rawSpikeBinsT,2),numNullDims);
    for nullDim = 1:numNullDims
        NullTraj(:,nullDim) = GPFATraj*NullSpace(:,nullDim);
    end
    
end