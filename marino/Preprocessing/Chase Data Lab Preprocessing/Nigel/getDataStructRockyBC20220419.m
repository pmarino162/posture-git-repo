function [Data] = getDataStructRockyBC20220419(rawData,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getSpikes = false;
        getSorts = false;
        exclCh = [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    getDecoder = true;
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = '';
    assignopts(who,varargin);
    
%% Unpack rawData
    taskStateCodes = rawData.TaskStateCodes;
    taskStateOutcomeMasks = rawData.TaskStateOutcomeMasks;
    outcomeMasks = rawData.OutcomeMasks;
    otherMasks = rawData.OtherMasks; %used for calib
    position = rawData.Position;
    velocity = rawData.Velocity;
    spikeCount = rawData.SpikeCount;
    time = rawData.Time;
    trialNo = rawData.TrialNo;
    
    rawSpikes = rawData.rawSpikes;

    
%% Preallocate Data Struct
    %stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    %Decoder
    Decoder = struct('decoderType','','runDecoder',[],'Parameters',[],'rawDecode',[],'rawSpikeBins',uint16([]),'name','','timestamps',[],'cursorTraj',[],'GPFATraj',[],'WorkTraj',[],'NullTraj',[]);
    %Data
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',[],'conditionData',[],'stateData',stateData,'Decoder',Decoder);
    
    %Marker
    if getMarker == true
        marker = struct('time',[],'position',[]);
        Data.marker = marker;
    end
    %Spikes
    if getSpikes == true
        %Load up rawSpikes Data for entire block
        channel = double(rawSpikes.channel);
        unit = double(rawSpikes.unit);
        source_index = double(rawSpikes.source_index);
        timestamps = double(rawSpikes.timestamps);
        binSourceTimestamps = double(rawSpikes.binSourceTimestamps);
        binSourceTimestampsAligned = double(rawSpikes.binSourceTimestampsAligned);       
        %Preallocated spikes struct
        spikes = struct('channel',[],'sort',[],'timestamps',[]);
        maxSort = max(unit);
        structInd = 1;
        for tempSource = 0:1
            for tempChannel = unique(channel)
                for tempSort = 0:1:maxSort
                    %If there are spikes for this channel-source-sort, preallocate
                    if any(source_index==tempSource & channel==tempChannel & unit==tempSort)
                           spikes(structInd).channel = tempChannel + 96*tempSource;
                           spikes(structInd).sort = tempSort;
                           spikes(structInd).timestamps = zeros(1,1000);
                           structInd = structInd + 1;
                    end
                end
            end
        end
        Data.spikes = spikes;
    end
    %Decoder
    if getDecoder == true
        Data.Decoder = struct('velocity',[],'position',[]);
    end
    numTrials = double(trialNo(end))+1;
    Data = repmat(Data,1,numTrials);

%% Fill Data Struct

%Fill Struct
for trial = 1:numTrials
    trial
    %Get trial inds from raw Data
        trialInds = find(trialNo==trial-1);
    %Get Trial Number
        Data(trial).trialNum = trial;
    %Get Trial Name
        Data(trial).trialName = trialName;
    %Get trial time - save this in decoder struct
        trialTime = time(trialInds)*1000;       %Convert to ms
        trialTime = trialTime - trialTime(1);    %Align to trial start
        
    %Get Trial Status
        Data(trial).trialStatus = outcomeMasks.Success(trialInds(1));
    %Get Target Data
        trialTargetData = position.target(1:2,trialInds);
        trialStateCodes = taskStateCodes.Values(trialInds);
        Data(trial).targetData = getTargetDataRocky(trialName,trialTargetData,trialStateCodes);
    %Get State Data
        Data(trial).stateData = getStateDataRocky(trialName,trialTime,trialStateCodes);
    %Get Spike Data
        if getSpikes == true
            trialBinSourceTimestamps = binSourceTimestamps(trialInds);
            trialStartTime = trialBinSourceTimestamps(1);
            trialEndTime = trialBinSourceTimestamps(end);
            for i = 1:numel(spikes)
                tempChannel = spikes(i).channel;
                if tempChannel > 96
                    searchChannel = tempChannel-96;
                    source = 1;
                else
                    searchChannel = tempChannel;
                    source = 0;
                end
                sort = spikes(i).sort;
                timestampsInd = channel==searchChannel & source_index==source & unit==sort;
                trialTimestamps = timestamps(timestampsInd);
                %Align, Convert to ms, round, and Save
                Data(trial).spikes(i).timestamps = round(1000.*(trialTimestamps(trialTimestamps>trialStartTime & trialTimestamps<=trialEndTime)-trialStartTime));
            end
        end  
    %Get Marker Data (Phasespace)
        if getMarker == true
        end
    %Get Decoder Data (if it's there)
        if getDecoder == true
            Data(trial).Decoder.spikes = spikeCount(:,trialInds)';
            Data(trial).Decoder.position = rawData.Position.Actual(1:2,trialInds)*1000;
            Data(trial).Decoder.velocity = rawData.Velocity.Decoded(1:2,trialInds)*1000;
            Data(trial).Decoder.time = round(trialTime);
        end
end

     
end