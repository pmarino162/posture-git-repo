function [Data] = getDataStructNigel20220419(rawData,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getSpikes = false;
        getSorts = false;
        exclCh = [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = '';
    assignopts(who,varargin);
    
%% Unpack rawData
    header = rawData.header;
    trials = rawData.trials;
    failTrials = rawData.fail_trials;
    catchTrials = rawData.catch_trials;
    rawSpikes = rawData.spikes;
    synch = rawData.synch;
    numTrials = size(synch.ComputerTime,1);
    feedback = rawData.em_feedback;

%% Reformat Chase Lab Spikes struct to match Batista Lab format
    spikes = struct('channel',[],'sort',[],'timestamps',[]);
    rawSpikesFieldNames = fieldnames(rawSpikes);
    spikes = repmat(spikes,1,size(rawSpikesFieldNames,1));
    for i = 1:size(rawSpikesFieldNames,1)
        %Get locations of channel and sort numbers within string
        firstChannelDigitInd = strfind(rawSpikesFieldNames{i,1},'Unit') + 4;
        lastChannelDigitInd = strfind(rawSpikesFieldNames{i,1},'_') - 1;
        firstSortDigitInd = strfind(rawSpikesFieldNames{i,1},'_') + 1;
        %Fill spikes struct with data
        spikes(i).timestamps = rawSpikes.(rawSpikesFieldNames{i})';
        spikes(i).channel = str2num(rawSpikesFieldNames{i,1}(firstChannelDigitInd:lastChannelDigitInd));
        spikes(i).sort = str2num(rawSpikesFieldNames{i,1}(firstSortDigitInd:end));
    end
    
%% Reformat trials structs 
    allTrials = trials([]);
    spikesStructShape = spikes([]);
    allTrials(1).spikes = spikesStructShape;
    allTrials = repmat(allTrials,1,numTrials);
    trialsFieldNames = fieldnames(trials);
    numFieldNames = size(trialsFieldNames,1);
    structInd = 1;
    %Successful trials
    for trialInd = 1:size(trials.Rewarded,1)
       for fieldInd = 1:numFieldNames
            allTrials(structInd).(trialsFieldNames{fieldInd})(1,:) = trials.(trialsFieldNames{fieldInd})(trialInd,:);
       end
       structInd = structInd + 1;
    end
    %Failed trials
    for trialInd = 1:size(failTrials.Rewarded,1)
       for fieldInd = 1:numFieldNames
            allTrials(structInd).(trialsFieldNames{fieldInd})(1,:) = failTrials.(trialsFieldNames{fieldInd})(trialInd,:);
       end
       structInd = structInd + 1;
    end
    %Catch trials
    for trialInd = 1:size(catchTrials.Rewarded,1)
       for fieldInd = 1:numFieldNames
            allTrials(structInd).(trialsFieldNames{fieldInd})(1,:) = catchTrials.(trialsFieldNames{fieldInd})(trialInd,:);
       end
       structInd = structInd + 1;
    end
    %Order trials chronologically 
    computerTrialTimes = [allTrials.ComputerTrialTime];    
    [~,sortInd] = sort(computerTrialTimes);
    allTrials = allTrials(sortInd);
        
%% Add each trial's spikes to allTrials
    %Trial duration measured by computer and plexon will be different. We
    %will use computer time to measure trial time, and we'll sync
    %computerTrialTime to plexonTrialTime throughout
    numSpikesRows = size(spikes,2);
    for trial = 1:numTrials
        %Get Trial Start and End Time
        computerStartTime = allTrials(trial).ComputerTrialTime;
        computerFinishTime = allTrials(trial).ComputerFinishTime;
        computerTrialDuration = computerFinishTime-computerStartTime;
        plexonStartTime = allTrials(trial).PlexonTrialTime;
        plexonFinishTime = allTrials(trial).PlexonFinishTime;
        %For each row of spikes, add relevant times to trial data
        for spikesRow = 1:numSpikesRows
            timestamps = spikes(spikesRow).timestamps;
            trialTimestampsMask = timestamps >= plexonStartTime & timestamps <= plexonStartTime + computerTrialDuration;
            timestamps = timestamps(trialTimestampsMask);
            channel = spikes(spikesRow).channel;
            sortNum = spikes(spikesRow).sort;
            allTrials(trial).spikes(spikesRow).timestamps = timestamps;
            allTrials(trial).spikes(spikesRow).channel = channel;
            allTrials(trial).spikes(spikesRow).sort = sortNum;
        end
    end

%% Add each trial's cursor positions to allTrials
    cursorTimestamps = feedback.Time;
    cursorPositions = feedback.CursorPosition;
    for trial = 1:numTrials
        %Get Trial Start and End Time
        computerStartTime = allTrials(trial).ComputerTrialTime;
        computerFinishTime = allTrials(trial).ComputerFinishTime;
        %For each row of spikes, add relevant times to trial data
        cursorTimestampsMask = (cursorTimestamps >= computerStartTime) & (cursorTimestamps <= computerFinishTime);
        allTrials(trial).cursor.timestamps = cursorTimestamps(cursorTimestampsMask);
        allTrials(trial).cursor.positions = cursorPositions(cursorTimestampsMask,1:2);
    end

%% Sync trials and adjust times; Convert to ms
    for trial = 1:numTrials
        %Convert trial events to ms
        allTrials(trial).HoldAStart = allTrials(trial).HoldAStart.*1000;
        allTrials(trial).HoldAFinish = allTrials(trial).HoldAFinish.*1000;
        allTrials(trial).HoldALength = allTrials(trial).HoldALength.*1000;
        allTrials(trial).ReactionFinish = allTrials(trial).ReactionFinish.*1000;
        allTrials(trial).HoldBStart = allTrials(trial).HoldBStart.*1000;
        allTrials(trial).HoldBFinish = allTrials(trial).HoldBFinish.*1000;
        allTrials(trial).HoldBLength = allTrials(trial).HoldBLength.*1000;
        %Sync Trial End Time, convert to ms
        computerStartTime = allTrials(trial).ComputerTrialTime;
        computerFinishTime = allTrials(trial).ComputerFinishTime;
        allTrials(trial).ComputerFinishTime = (computerFinishTime-computerStartTime).*1000;
        %Sync Spikes, convert to ms
        plexonTrialTime = allTrials(trial).PlexonTrialTime;
        spikes = allTrials(trial).spikes;
        for spikesRow = 1:numSpikesRows
            timestamps = spikes(spikesRow).timestamps;
            spikes(spikesRow).timestamps = (timestamps-plexonTrialTime).*1000;
        end
        allTrials(trial).spikes = spikes;
        %Sync cursor positions, convert to ms
        allTrials(trial).cursor.timestamps = (allTrials(trial).cursor.timestamps - computerStartTime).*1000;
    end

%% Convert Lengths to mm
    %Header
        header.cout.CursorRadius = header.cout.CursorRadius*1000;
        header.cout.TargetRadius = header.cout.TargetRadius*1000;
    %Trials
        for trial = 1:numTrials
            allTrials(trial).StartPos = allTrials(trial).StartPos.*1000;
            allTrials(trial).TargetPos = allTrials(trial).TargetPos.*1000;
            allTrials(trial).cursor.positions = allTrials(trial).cursor.positions.*1000;
        end
        
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
        Data.spikes = struct('channel',[],'sort',[],'timestamps',[]);
    end

    Data = repmat(Data,1,numTrials);

%% Fill Data Struct
for trial = 1:numTrials
    %Get Trial Number
        Data(trial).trialNum = trial;
    %Get Trial Name
        Data(trial).trialName = trialName;
    %Get Trial Status
        Data(trial).trialStatus = allTrials(trial).Rewarded;
    %Get Target Data
        Data(trial).targetData = getTargetDataNigel(trialName,header,allTrials(trial));
    %Get State Data (and convert transitions to ms)
        Data(trial).stateData = getStateDataNigel(trialName,allTrials(trial));
    %Get Spike Data
        if getSpikes == true
            spikes = allTrials(trial).spikes;
            if ~isempty(spikes)
                %Exlcude channels
                spikes(ismember([spikes.channel],exclCh)) = [];
                %If ~getSorts, combine across channels
                if ~getSorts
                    tempSpikes = struct('channel',[],'sort',[],'timestamps',[]);
                    uniqueChannels = unique([spikes.channel]);
                    tempSpikesInd = 1;
                    for channel = uniqueChannels
                       tempSpikes(tempSpikesInd).channel = channel;
                       tempSpikes(tempSpikesInd).sort = 'all';
                       timestamps = sort([spikes([spikes.channel]==channel).timestamps]);
                       tempSpikes(tempSpikesInd).timestamps = timestamps;
                       tempSpikesInd = tempSpikesInd + 1; 
                    end
                    spikes  = tempSpikes;
                end
            end
            %Store
            Data(trial).spikes = spikes;
        end
    %Get Marker Data (Phasespace)
        if getMarker == true
            if ~isempty(rawPositions)
                %Center Marker Positions
            end
        end

    %Get Decoder Data (if it's there)
        %Interpolate Cursor Position up to 1000Hz 
        time = allTrials(trial).cursor.timestamps;
        cursorTraj = allTrials(trial).cursor.positions;
        interpTime = ceil(time(1)):1:floor(time(end));
        cursorTraj = interp1(time,cursorTraj,interpTime);
        Data(trial).Decoder.timestamps = interpTime;
        Data(trial).Decoder.cursorTraj = cursorTraj;
end
     
end