function [Data] = getDataStruct(rawData,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getMarker = false;
    getKin = false;
    checkPhasespaceSync = false;
    getForce = false;
    getAlg = false;
    getEye = false;
    getSpikes = false;
    getSorts = false;
    getSmoothedFR = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    inclStateTable = false;
    removeBadPhasespace = true;
    binWidth = 1;
    exclCh = [];
    exclTrials = [];
    exclZero = true; %Exclude channels that only contain sort 0
    centerMarker = true; %Subracts workspace center from marker pose if true
    forceSetup = '';
    assignopts(who,varargin);
    
%% Preprocess Data
    procData = preprocessDataMarino(rawData,'CHECKPHASESPACESYNC',checkPhasespaceSync,'DROPPEDMARKERTRIALS','all','CHECKDECODE',false);
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
    
%% Preallocate Data Struct
    numTrials = size(procData,2);
    %stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    %Decoder
    Decoder = struct('decoderType','','runDecoder',[],'Parameters',[],'rawDecode',[],'rawSpikeBins',uint16([]),'name','','timestamps',[],'cursorTraj',[],'GPFATraj',[],'WorkTraj',[],'NullTraj',[]);
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
        leftPupilInd = find(strcmpi(analogChannelNames,'Left Pupil'));
        rightPupilInd = find(strcmpi(analogChannelNames,'Right Pupil'));
    end
    %Spikes
    if getSpikes == true
%         spikes = procData(1).TrialData.spikes;
%         [spikes] = preallocateSpikesStruct(spikes,exclZero,getSorts)
        
        channelSpikes = struct('channel',uint16([]),'sort',uint16([]),'binnedSpikes',uint16([]));
        spikes = procData(1).TrialData.spikes;
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
        Data.spikes = spikes;
      
%         
%         if getSmoothedFR == true
%             channelSpikes = struct('channel',uint16([]),'sort',uint16([]),'binnedSpikes',uint16([]),'smoothedFR',double([]));
%         else
%             channelSpikes = struct('channel',uint16([]),'sort',uint16([]),'binnedSpikes',uint16([]));
%         end
        
%         if getSmoothedFR == true
%             channelSpikes = struct('sorts',struct('sortID',uint16([]),'bins',uint16([]),'smoothedFR',double([])),'allSorts',uint16([]));
%         else
%             channelSpikes = struct('sorts',struct('sortID',uint16([]),'bins',uint16([])),'allSorts',uint16([]));
%         end
%         channelSpikes = repmat(channelSpikes,1,numChannels);
%         for channel = 1:numChannels
%             chSpikes = spikes([spikes.channel]==channel);
%             numSorts(channel) = size(chSpikes,2);   
%             for sortInd = 1:numSorts(channel)
%                 channelSpikes(channel).sorts(sortInd).sortID = chSpikes(sortInd).sort;
%             end
% %         end
%         if getSmoothedFR == true
%             spikes = struct('binTimes',uint32([]),'channelSpikes',channelSpikes,'allChannelSpikeBins',uint16([]),'allSortSpikeBins',uint16([]),'allChannelSmoothedFR',double([]),'allSortSmoothedFR',double([]));
%         else
%             spikes = struct('binTimes',uint32([]),'channelSpikes',channelSpikes,'allChannelSpikeBins',uint16([]),'allSortSpikeBins',uint16([]));
%         end
%         Data.spikes = spikes;
    end

    Data = repmat(Data,1,numTrials);
    
%% Fill Data Struct
    for trial = 1:numTrials
%         trial
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
            stateTransitions = procData(trial).TrialData.stateTransitions;
            time = getActualTime(procData(trial).Definitions.analogChannelNames,procData(trial).TrialData.analogData);
            %Convert transition frames to transition times
            for i = 1:size(stateTransitions,2)
               stateTransitions(2,i) = time(stateTransitions(2,i),1);
            end
            Data(trial).stateData(1).stateNames = stateNames;
            Data(trial).stateData(1).stateTransitions = stateTransitions;
            if inclStateTable == true
               Data(trial).stateData(1).stateTable = procData(trial).Parameters.StateTable; 
            end
        %Get Marker Data (Phasespace)
            if getMarker == true
                rawPositions = procData(trial).TrialData.Marker.rawPositions;
                if ~isempty(rawPositions)
                    time = rawPositions(:,6)';
                    position = rawPositions(:,2:3);
                    %Flip Y
                    position(:,2) = -1.*position(:,2);
                    %Center Marker Positions
                    if centerMarker
                        workspaceCenter = Data(trial).targetData.workspaceCenter(1,1:2);
                        position = position - workspaceCenter;
                    end
                    %Filter position
                    position = filtfilt(b,a,position);
                    %Differentiate position to get velocity
                    dp = diff(position);
                    dt = diff(time);
                    velocity = dp./dt';
                    velocityTime = time(1:end-1)+dt/2;
                    %Interpolate both signals up to 1000Hz
%                     interpTime = ceil(time(1)):1:floor(time(end));
%                     position = interp1(time,position,interpTime);
%                     velocity = interp1(velocityTime,velocity,interpTime);
                    %Add to data struct
                    Data(trial).marker.time = time;
                    Data(trial).marker.position = position;
                    Data(trial).marker.velocity = velocity;
                end
            end
        %Get Force Data
            if getForce == true
                rawForces = procData(trial).TrialData.Marker.rawForces;
                forceCursor = procData(trial).TrialData.Marker.forceCursor;
                allFTOutput = procData(trial).TrialData.Marker.allFTOutput;
                absoluteForcesAndTorques = procData(trial).TrialData.Marker.absoluteForcesAndTorques;
%                 %Interpolate all signals up to 1000Hz
%                 time = rawForces(:,1);
%                 interpTime = ceil(time(1)):1:floor(time(end));
%                 force = interp1(time,rawForces(:,2:5),interpTime);
%                 forceCursor = interp1(time,forceCursor(:,2:4),interpTime);
%                 allFTOutput = interp1(time,allFTOutput,interpTime);
%                 absoluteForcesAndTorques = interp1(time,absoluteForcesAndTorques(:,2:7),interpTime);
               %
                Data(trial).force.time = rawForces(:,1)';
                Data(trial).force.force = rawForces(:,2:5);
                Data(trial).force.forceCursor = forceCursor(:,2:4);
                Data(trial).force.allFTOutput = allFTOutput;
                Data(trial).force.absoluteForcesAndTorques = absoluteForcesAndTorques(:,2:7);
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
                    pupil = [analogData(:,leftPupilInd),analogData(:,rightPupilInd)];
                    %Flip Y
                    position(:,2) = -1.*position(:,2);
                    %Center Eye Positions
                    workspaceCenter = Data(trial).targetData.workspaceCenter(1,1:2);
                    position = position - workspaceCenter;
                    %Store Data
                    Data(trial).eye.time = time;
                    Data(trial).eye.position = position;
                    Data(trial).eye.pupil = pupil;
                end
            end
        %Get Spike Data (optional)
            if getSpikes == true
                channelSpikes = Data(trial).spikes.channelSpikes;
                spikes = procData(trial).TrialData.spikes;
                if getSorts == true
                    [binTimes,channelSpikes,~,allSortSpikeBins] = getBinnedSpikes(binWidth,spikes,stateTransitions,channelSpikes,getSorts);
                    Data(trial).spikes.allSortSpikeBins = allSortSpikeBins;
                else
                    [binTimes,channelSpikes,allChannelSpikeBins,~] = getBinnedSpikes(binWidth,spikes,stateTransitions,channelSpikes,getSorts);
                    Data(trial).spikes.allChannelSpikeBins = allChannelSpikeBins;
                end
                Data(trial).spikes.binTimes = binTimes;
                Data(trial).spikes.channelSpikes = channelSpikes;
%                 [binTimes,channelSpikes,allChannelSpikeBins,allSortSpikeBins] = getBinnedSpikes(binWidth,spikes,stateTransitions,channelSpikes,exclCh);
%                 %Exclude Channels
%                     %allSortSpikeBins
%                     allSortSpikeBinsColumnInd = 1;
%                     allSortSpikeBinsColumnsToDelete = [];
%                     for channel = 1:numChannels
%                        numChSorts = numSorts(channel);
%                        if ismember(channel,exclCh)
%                             allSortSpikeBinsColumnsToDelete  = [allSortSpikeBinsColumnsToDelete,allSortSpikeBinsColumnInd:allSortSpikeBinsColumnInd+numChSorts-1];
%                        end
%                        allSortSpikeBinsColumnInd = allSortSpikeBinsColumnInd + numChSorts;
%                     end
%                     allSortSpikeBins(:,allSortSpikeBinsColumnsToDelete) = [];
%                     %channelSpikes
%                     channelSpikes(exclCh) = [];
%                     %allChannelSpikeBins
%                     allChannelSpikeBins(:,exclCh) = [];
                %Store Data
%                 Data(trial).spikes.binTimes = binTimes;
%                 Data(trial).spikes.channelSpikes = channelSpikes;
%                 Data(trial).spikes.allChannelSpikeBins = allChannelSpikeBins;
%                 Data(trial).spikes.allSortSpikeBins = allSortSpikeBins;

%                 %Get Smoothed FR (optional) - This functionality has been
%                 replaced/improved in getStatesTraj
%                 if getSmoothedFR == true
%                    %allChannelSpikeBins
%                    allChannelSpikeBins = double(Data(trial).spikes.allChannelSpikeBins); 
%                    allChannelSmoothedFR = allChannelSpikeBins./binWidth.*1000;
%                    allChannelSmoothedFR = smoothdata(allChannelSmoothedFR,1,'gaussian',5);
%                    Data(trial).spikes.allChannelSmoothedFR = allChannelSmoothedFR;
%                    %allSortSpikeBins
%                    allSortSpikeBins = double(Data(trial).spikes.allSortSpikeBins); 
%                    allSortSmoothedFR = allSortSpikeBins./binWidth.*1000;
%                    allSortSmoothedFR = smoothdata(allSortSmoothedFR,1,'gaussian',5);
%                    Data(trial).spikes.allSortSmoothedFR = allSortSmoothedFR;
%                 end
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
%                 rawDecode = Decoder.rawDecode;
                Data(trial).Decoder.name = Decoder.Parameters.name;
                Data(trial).Decoder.timestamps = rawDecode(:,1)';
                Data(trial).Decoder.cursorTraj = rawDecode(:,3:4);
                if strcmpi(decoderType,'GPFA')
                    [GPFATraj,WorkTraj,NullTraj] = getDecoderGPFATraj(Decoder);
                    Data(trial).Decoder.GPFATraj = GPFATraj;
                    Data(trial).Decoder.WorkTraj = WorkTraj;
                    Data(trial).Decoder.NullTraj = NullTraj;
                    %Flip Y
                    Data(trial).Decoder.WorkTraj(:,2) = -Data(trial).Decoder.WorkTraj(:,2);
                end
                %Flip Y
                Data(trial).Decoder.cursorTraj(:,2) = -Data(trial).Decoder.cursorTraj(:,2);
            end
            %Get Kinematic Data
            if getKin == true
                Data(trial).kinData = getKinData(trialName,procData(trial),Data(trial),centerMarker);
            end
    end
    
    %Exclude Trials
    trialNum = [Data.trialNum];
    rmList = find(ismember(trialNum,exclTrials));
    Data(rmList) = [];
    
end