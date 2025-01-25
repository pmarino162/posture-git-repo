function [Data] = getDataStructRockyBC_EZ_20220419(EZ,varargin)

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
    

%% Preallocate Data Struct
    %stateData
    stateData = struct('stateNames',{},'stateTransitions',[]);
    %targetData
    targetData = struct('centerLoc',[],'centerSize',[],'targetLoc',[],'targetSize',[],'targetID',[]);
    %Decoder
    Decoder = struct('decoderType','','runDecoder',[],'Parameters',[],'rawDecode',[],'rawSpikeBins',uint16([]),'name','','timestamps',[],'cursorTraj',[],'GPFATraj',[],'WorkTraj',[],'NullTraj',[]);
    %Data
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',targetData,'conditionData',[],'stateData',stateData,'Decoder',Decoder);
    %Marker
    if getMarker == true
        marker = struct('time',[],'position',[]);
        Data.marker = marker;
    end
    %Spikes
    if getSpikes == true
        %Load up rawSpikes Data for entire block
        spike_data_1 = EZ.spike_data_1;
        numChannels = size(spike_data_1,1);
        numSorts = size(spike_data_1,2)-1;
        %Preallocated spikes struct
        spikes = struct('channel',[],'sort',[],'timestamps',[]);
        structInd = 1;
        for channel = 1:numChannels
            for sort = 1:numSorts
                spikes(structInd).channel = channel;
                spikes(structInd).sort = sort;
                spikes(structInd).timestamps = zeros(1,1000);
                structInd = structInd + 1;
            end
        end
        Data.spikes = spikes;
    end
    %Decoder
    if getDecoder == true
        Data.Decoder = struct('velocity',[],'position',[]);
    end
    %Get number of trials; preallocate Data
    fields = fieldnames(EZ);
    numTrials = str2num(erase(fields{end},'trial_header_'));
    Data = repmat(Data,1,numTrials);

%% Get List of all possible state names; assign them numbers
    %Use commented code to see all states, then manually create list, since
    %states could have appeared in a weird order
%     allStateNames = {}; allStateNamesInd = 1;
%     for trial = 1:numTrials
%         trialHeaderName = ['trial_header_',num2str(trial)];
%         trial_header = EZ.(trialHeaderName);
%         trialStateNames = trial_header.beh_event.type;
%         for i = 1:numel(trialStateNames)
%             state = trialStateNames{1,i};
%             if ~any(cellfun(@(x) strcmp(state,x), allStateNames))
%                 allStateNames{1,allStateNamesInd} = state;
%                 allStateNamesInd = allStateNamesInd + 1;
%             end
%         end
%     end
    stateNames = {'SessionStart','Center','HoldA','React','Move','Hold','InterTrial','TrialEnd'};
    

%% Fill Data Struct
%Fill Struct
for trial = 1:numTrials
    %Unpack trial info
        trialHeaderName = ['trial_header_',num2str(trial)];
        trial_header = EZ.(trialHeaderName);
        spikeDataName = ['spike_data_',num2str(trial)];
        spike_data = EZ.(spikeDataName);
        extractionHeaderName = ['extraction_header_',num2str(trial)];
        extraction_header = EZ.(extractionHeaderName);
        extractionDataName = ['extraction_data_',num2str(trial)];
        extraction_data = EZ.(extractionDataName);
        
    %Get Trial Number
        Data(trial).trialNum = trial;
    %Get Trial Name
        Data(trial).trialName = trialName;       
    %Get Trial Status
        Data(trial).trialStatus = trial_header.success;
    %Get Target Data
        targetData.centerLoc = [0,0];
        targetData.centerSize = trial_header.target.radii.center*1000;
        targetData.targetSize = trial_header.target.radii.periph*1000;
        targetID = trial_header.target_id;
            if targetID==5
                targetID = 2;
            elseif targetID == 6
                targetID = 4;
            elseif targetID == 2
                targetID = 5;
            elseif targetID ==7
                targetID = 6;
            elseif targetID == 4
                targetID = 7;
            end
        targetData.targetID = targetID;
        targetData.targetLoc = 1000*trial_header.target.distance*[cosd((trial_header.target_id-1)*45), sind((trial_header.target_id-1)*45)];
        Data(trial).targetData = targetData;
     %Get State Data
        stateData(1).stateNames = stateNames;
        trialStateNames = trial_header.beh_event.type;
        numTrialStates = numel(trialStateNames);
        stateTransitions = zeros(2,numTrialStates);
        for trialState = 1:numTrialStates
           stateTransitions(1,trialState) = find(cellfun(@(x) strcmp(trialStateNames{1,trialState},x), stateNames)); 
        end
        %Align to state transitions to start, convert to ms, round
        stateTransitions(2,:) = round(1000*(trial_header.beh_event.time-trial_header.beh_event.time(1,1)))';
        stateData(1).stateTransitions = stateTransitions;
        Data(trial).stateData = stateData;
     %Get Spike Data
        if getSpikes == true
            for i = 1:numel(spikes)
               channel = spikes(i).channel;
               sort = spikes(i).sort;
               %Align to trial start; convert to ms
               Data(trial).spikes(i).timestamps = (1000*(spike_data{channel,sort+1}-trial_header.beh_event.time(1,1)))';
            end
        end  
    %Get Marker Data (Phasespace)
        if getMarker == true
        end
    %Get Decoder Data (if it's there)
        if getDecoder == true
            pos = extraction_data.pos';
            vel = extraction_data.vel';
            sernum_vel = extraction_header.sernum_vel;
            sernum_pos = extraction_header.sernum_pos;
            send_time_vel = extraction_header.send_time_vel;
            send_time_pos = extraction_header.send_time_pos;
            dt = extraction_header.dt;
            trialStartTimeComputer = trial_header.task_state.time(1);
            %Sync and convert to ms
            posTime = round(1000*(send_time_pos-trialStartTimeComputer));
            velTime = round(1000*(send_time_vel-trialStartTimeComputer));
            %Remove any repeated samples
            [posTime,uniquePosInd,~] = unique(posTime,'first');
            pos = pos(uniquePosInd,:);         
            [velTime,uniqueVelInd,~] = unique(velTime,'first');
            vel = vel(uniqueVelInd,:);
            %Interpolate decoder data up to 1000Hz
            interpPosTime = ceil(posTime(1)):1:floor(posTime(end));
            interpVelTime = ceil(velTime(1)):1:floor(velTime(end));
            pos = interp1(posTime,pos,interpPosTime);
            vel = interp1(velTime,vel,interpVelTime);   
            %Store
            Data(trial).Decoder.spikes = [];
            Data(trial).Decoder.position = pos;
            Data(trial).Decoder.velocity = vel;
            Data(trial).Decoder.posTime = interpPosTime;
            Data(trial).Decoder.velTime = interpVelTime;
                %send_time = Int.QL.Headers.TASK_STATE_CONFIG.send_time; %One of these for every state
                %Code that I used to test that there was one send time for
                %every state and not one for every "behavioral event"
                %                 fields = fieldnames(EZ);
                %                 numStates = 0;
                %                 numBE = 0;
                %                 for i = 1:size(fields,1)
                %                    trial = i;
                %                    trialHeaderName = ['trial_header_',num2str(trial)];
                %                    trial_header = EZ.(trialHeaderName);
                %                    numTrialStates = size(trial_header.task_state.id,2);
                %                    numStates = numStates + numTrialStates;
                %                    numTrialBE = size(trial_header.beh_event.type,2);
                %                    numBE = numBE + numTrialBE;
                %                 end
            

        end
end

     
end