function [Data] = getDataStruct20211210(rawData,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getSpikes = false;
        getSorts = false;
        exclCh = [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = true;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    getForce = false;
        forceSetup = '';
        centerForceCursor = true; %Subracts workspace center from force cursor pose if true
    getAlg = false;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = false;
    checkDecode = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];
    assignopts(who,varargin);
    
%% Preprocess Data
    procData = preprocessDataMarino(rawData,'CHECKPHASESPACESYNC',checkPhasespaceSync,'DROPPEDMARKERTRIALS',droppedMarkerTrials,'CHECKDECODE',checkDecode,'REMOVESCREENFREEZETRIALS',removeScreenFreezeTrials);
    clearvars rawData
    
    %if there are duplicate psp samples, remove them
    numTrials = size(procData,2);
    for trial = 1:numTrials
        rowsToDelete = [];
        rawPositions = procData(trial).TrialData.Marker.rawPositions;
        if ~isempty(rawPositions)
            time = rawPositions(:,6);
            numRows = size(rawPositions,1);
            for row = 1:numRows-1
                if time(row+1,1) == time(row,1)
                    rowsToDelete = [rowsToDelete,row+1];
                end
            end
            rawPositions(rowsToDelete,:) = [];
            procData(trial).TrialData.Marker.rawPositions = rawPositions;
        end
    end

    %remove any trials for which the phasespace data wasn't saved properly
    if removeBadPhasespace
        numTrials = size(procData,2);
        rmTrials = [];
        for trial = 1:numTrials
            rawPositions = procData(trial).TrialData.Marker.rawPositions;
            if ~isempty(rawPositions)
                pspEndTime = rawPositions(end,6);
                trialEndTime = procData(trial).TrialData.stateTransitions(2,end);
                if trialEndTime > pspEndTime + 500
                    rmTrials = [rmTrials,trial];
                end
            end
        end
        procData(rmTrials) = [];
    end
    
%% Remove Catch Trials
    procData = removeCatchTrials(procData);
    
%% Remove Exclude Trials
    trialNumber = NaN(1,numel(procData));
    for i = 1:numel(procData)
        curTrialNumber = procData(i).Overview.trialNumber;
        trialNumber(i) = str2num(curTrialNumber(6:end));
    end
    rmList = find(ismember(trialNumber,exclTrials));
    procData(rmList) = [];
    
%% Preallocate Data Struct
    numTrials = size(procData,2);
    %stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    %Decoder
    %Decoder = struct('decoderType','','runDecoder',[],'Parameters',[],'rawDecode',[],'rawSpikeBins',uint16([]),'name','','timestamps',[],'cursorTraj',[],'GPFATraj',[],'WorkTraj',[],'NullTraj',[]);
    Decoder = struct('decoderType','','Parameters',[],'rawDecode',[],'rawSpikeBins',uint16([]),'name','','timestamps',[],'cursorTraj',[]);
    %Data
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',[],'conditionData',[],'stateData',stateData,'Decoder',Decoder);

    %Marker
    if getMarker == true
        marker = struct('time',[],'position',[],'velocity',[]);
        Data.marker = marker;
        %Design butterworth filter for phasepace data
        rawPositions = procData(1).TrialData.Marker.rawPositions;
        time = rawPositions(:,6)';
        fs =  1/(mode(diff(time))/1000);                       %sampling frequency (Hz)
        Wn = 15/fs;                                            %cutoff frequency (normalized)
        [b,a] = butter(4,Wn);                                  %4th order butterworth
    end
    %Force
    if getForce == true
        force = struct('time',[],'force',[],'torque',[],'forceCursor',[]);
        Data.force = force;
        procData = getForceAndForceCursor(procData,forceSetup);
    end
    %Analog Data 
    if getAlg == true
        Data.algData = [];
    end
    %Eye
    if getEye == true
        eye = struct('time',[],'position',[],'pupil',[]);
        Data.eye = eye;
        procData = calibrateEyePosition(procData);
        analogChannelNames = procData(1).Definitions.analogChannelNames;
        leftPupilInd = find(strcmpi(analogChannelNames,'Left Pupil'));
        rightPupilInd = find(strcmpi(analogChannelNames,'Right Pupil'));
    end
    %Spikes
    if getSpikes == true
        Data.spikes = struct('channel',[],'sort',[],'timestamps',[]);
        %Identify zero ch for exlcusion (optional)
        if exclZero
            zeroCh = [];
            spikes = procData(1).TrialData.spikes;
            for channel = unique([spikes.channel])                        
                chSorts = [spikes([spikes.channel]==channel).sort];
                if size(chSorts,2) == 1 && chSorts(1) == 0
                    zeroCh = [zeroCh,channel];
                end
            end
        end
        zeroCh = setdiff(zeroCh,exclCh);
    end

    Data = repmat(Data,1,numTrials);
    
%% Fill Data Struct
    for trial = 1:numTrials
        trial
        %Get Trial Number
            trialNumber = procData(trial).Overview.trialNumber;
            trialNumber = trialNumber(6:end);
            Data(trial).trialNum = str2double(trialNumber);
        %Get Trial Name
            trialName = procData(trial).Overview.trialName;
            Data(trial).trialName = trialName;
        %Get Trial Status
            Data(trial).trialStatus = procData(trial).Overview.trialStatus;
        %Get Target Data (if it's there)
            if ~isempty(procData(trial).Parameters.TrialTargets)
                TrialTargets = procData(trial).Parameters.TrialTargets;
                Data(trial).targetData = getTargetData(trialName,TrialTargets,procData(trial));
            end
        %Get Condition Data
            Data(trial).conditionData = getConditionData(trialName,procData(trial),Data(trial));
        %Get Tube Data (if they were used)
            tubeWindowCheckType = double([procData(trial).Parameters.StateTable.tubeWindowCheckType]);
            if sum(tubeWindowCheckType ~= 0) > 0
                TrialTubeParameters = procData(trial).Parameters.TrialTubeParameters;
                Data(trial).tubeData.tubeName = TrialTubeParameters.tubeName;
                Data(trial).tubeData.trajectory = TrialTubeParameters.trajectory;
                %Flip Y
                Data(trial).tubeData.trajectory(:,2) = -1.*Data(trial).tubeData.trajectory(:,2);
                window = TrialTubeParameters.window(:,1)';
                windowMask = isnan(window)==0;
                Data(trial).tubeData.size = unique(window(windowMask));
            end
        %Get State Data
            stateNames = procData(trial).Parameters.stateNames;
            stateTransitions = double(procData(trial).TrialData.stateTransitions);
            time = getActualTime(procData(trial).Definitions.analogChannelNames,procData(trial).TrialData.analogData);
            %Convert transition frames to transition times
            for i = 1:size(stateTransitions,2)
               stateTransitions(2,i) = time(stateTransitions(2,i),1);
            end
            Data(trial).stateData(1).stateNames = stateNames;
            Data(trial).stateData(1).stateTransitions = round(stateTransitions);
            if inclStateTable == true
               Data(trial).stateData(1).stateTable = procData(trial).Parameters.StateTable; 
            end
        %Get Marker Data (Phasespace)
            if getMarker == true
                rawPositions = procData(trial).TrialData.Marker.rawPositions;
                if ~isempty(rawPositions)
                    time = rawPositions(:,6)';
                    position = rawPositions(:,2:4);
                    %Flip Y
                    position(:,2) = -1.*position(:,2);
                    %Center Marker Positions
                    if centerMarker
                        workspaceCenter = Data(trial).targetData.workspaceCenter;%(1,1:2);
                        position = position - workspaceCenter;
                    end
                    %Filter position
                    position = filtfilt(b,a,position);
                    %Differentiate position to get velocity (mm/ms)
                    dp = diff(position);
                    dt = diff(time);
                    velocity = dp./dt';
                    velocityTime = time(1:end-1)+dt/2;
                    %Interpolate both signals up to 1000Hz
                    interpTime = ceil(time(1)):1:floor(time(end));
                    position = interp1(time,position,interpTime);
                    velocity = interp1(velocityTime,velocity,interpTime);
                    %Add to data struct
                    Data(trial).marker.time = interpTime;
                    Data(trial).marker.position = position;
                    Data(trial).marker.velocity = velocity;
                end
            end
        %Get Force Data
            if getForce == true
                %Get data
                rawForces = procData(trial).TrialData.Marker.rawForces;
                forceCursor = procData(trial).TrialData.Marker.forceCursor;
                allFTOutput = procData(trial).TrialData.Marker.allFTOutput;
                absoluteForcesAndTorques = procData(trial).TrialData.Marker.absoluteForcesAndTorques;
                %Flip Force Cursor Y
                forceCursor(:,3) = -1.*forceCursor(:,3);
                %Center Marker Positions
                if centerForceCursor
                    workspaceCenter = Data(trial).targetData.workspaceCenter(1,1:2);
                    forceCursor(:,2:3) = forceCursor(:,2:3) - workspaceCenter;
                end
                %Interpolate all signals up to 1000Hz
                time = rawForces(:,1);
                interpTime = ceil(time(1)):1:floor(time(end));
                force = interp1(time,rawForces(:,2:5),interpTime);
                forceCursor = interp1(time,forceCursor(:,2:4),interpTime);
                allFTOutput = interp1(time,allFTOutput,interpTime);
                absoluteForcesAndTorques = interp1(time,absoluteForcesAndTorques(:,2:7),interpTime);
                %Store
                Data(trial).force.time = interpTime;
                Data(trial).force.force = force;
                Data(trial).force.forceCursor = forceCursor;
                Data(trial).force.allFTOutput = allFTOutput;
                Data(trial).force.absoluteForcesAndTorques = absoluteForcesAndTorques;
            end
        %Get Analog Data
            if getAlg == true
               algData = procData(trial).TrialData.analogData;
               Data(trial).algData = algData;
            end
        %Get Eye Data
            if getEye == true
                position = procData(trial).TrialData.EyeData.position;
                if size(position,2) > 1 
                    time = getActualTime(procData(trial).Definitions.analogChannelNames,procData(trial).TrialData.analogData)';
                    position = position(:,1:2);
                    pupil = [algData(:,leftPupilInd),algData(:,rightPupilInd)];
                    %Flip Y
                    position(:,2) = -1.*position(:,2);
                    %Center Eye Positions
                    workspaceCenter = Data(trial).targetData.workspaceCenter(1,1:2);
                    position = position - workspaceCenter;
                    %Interpolate to 1000Hz
                    interpTime = ceil(time(1)):1:floor(time(end));
                    position = interp1(time,position,interpTime);
                    pupil = interp1(time,pupil,interpTime);
                    %Store Data
                    Data(trial).eye.time = interpTime;
                    Data(trial).eye.position = position;
                    Data(trial).eye.pupil = pupil;
                end
            end
        %Get Spike Data (optional)
            if getSpikes == true
                spikes = procData(trial).TrialData.spikes;
                %Convert timestamps field to double
                for i = 1:size(spikes,2)
                    spikes(i).timestamps = double(spikes(i).timestamps);
                end
                if ~isempty(spikes)
                    %Delete source and timestampUnits fields; reorder fields
                    spikes = rmfield(spikes,{'source','timestampUnits'});
                    if isfield(spikes,'spikeWaveforms')
                        spikes = rmfield(spikes,{'spikeWaveforms'});
                    end
                    spikes = orderfields(spikes,Data(trial).spikes);
                    %Exlcude channels
                    spikes(ismember([spikes.channel],exclCh)) = [];
                    %Exclude channels that only contain zero sorts (optional)
                    if exclZero
                        spikes(ismember([spikes.channel],zeroCh)) = [];
                    end
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
        %Get Waveform Data (optional)
            if getWaveforms == true
                Data(trial).waveforms = getWaveformData(spikes);
            end
        %Get LFP Data (optional)
            if getLFP == true
               lfp = procData(trial).TrialData.lfp;
               Data(trial).LFP = lfp;
            end
        %Get Broadband Data (optional)
            if getBB == true
               bb_waveforms = double(procData(trial).TrialData.TDT.bb_waveforms);
               bb_samplingRate = double(procData(trial).TrialData.TDT.bb_samplingRate);
               BB = struct('waveforms',bb_waveforms,'fs',bb_samplingRate);
               Data(trial).BB = BB;
            end
        %Get Decoder Data (if it's there)
            Decoder = procData(trial).TrialData.Decoder;
            decoderType = Decoder.decoderType;
            rawDecode = Decoder.rawDecode;
            if ~strcmpi(decoderType,'NULL') && ~isempty(rawDecode)
                Data(trial).Decoder.decoderType = Decoder.decoderType;
                Data(trial).Decoder.runDecoder = Decoder.runDecoder;
                Data(trial).Decoder.Parameters = Decoder.Parameters;
                Data(trial).Decoder.rawSpikeBins = uint16(Decoder.rawSpikeBins);
                %Interpolate Cursor Position up to 1000Hz 
                time = rawDecode(:,1)';
                cursorTraj = rawDecode(:,3:4);
                interpTime = ceil(time(1)):1:floor(time(end));
                cursorTraj = interp1(time,cursorTraj,interpTime);
                Data(trial).Decoder.name = Decoder.Parameters.name;
                Data(trial).Decoder.timestamps = interpTime;
                Data(trial).Decoder.cursorTraj = cursorTraj;
                Data(trial).Decoder.decoderTime = rawDecode(:,1)';
                Data(trial).Decoder.decoderTraj = rawDecode(:,3:4);
                %if strcmpi(decoderType,'GPFA')
                %    [GPFATraj,WorkTraj,NullTraj] = getDecoderGPFATraj(Decoder);
                %    Data(trial).Decoder.GPFATraj = GPFATraj;
                %    Data(trial).Decoder.WorkTraj = WorkTraj;
                %    Data(trial).Decoder.NullTraj = NullTraj;
                %    %Flip Y
                %    Data(trial).Decoder.WorkTraj(:,2) = -Data(trial).Decoder.WorkTraj(:,2);
                %end
                %Flip Y
                Data(trial).Decoder.cursorTraj(:,2) = -Data(trial).Decoder.cursorTraj(:,2);
                Data(trial).Decoder.decoderTraj(:,2) = -Data(trial).Decoder.decoderTraj(:,2);
            end
        %Get Kinematic Data - REPLACED W FXN AT END OF
        %preprocessAndSaveData20220419
            %if getKin == true
                %Data(trial).kinData = getKinData(trialName,procData(trial),Data(trial),centerMarker);
            %end
    end
    

    
end