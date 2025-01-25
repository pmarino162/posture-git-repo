function [TD] = getTrajectoryData(T,stateNames,targetState,varargin)
% Return trajectory and neural data for specified state
%
% [T] = getTrajectoryData(T,stateNames,targetState,varargin)
%
% Inputs (required):
%   T           Preprocessed data structure (Trial format)
%   stateNames  Cell array of states to include
%   targetState State containing target location.
%
% Outputs:
%   TD          Trajectory data object

% Optional arguments:
% Align to movement onset?  Currently, t=0 corresponds to the onset of the
% trajectory state.  However, in cases where movement onset is defined, it
% would be desirable to align to movement onset.

% Set default parameters for optional arguments
alignMoveOnset      = false;    % Specify whether to align to movment onset (not currently used)
successOnly         = true;     % Remove invalid trails
rmNoSpikeData       = true;     % Remove trials with no spike data
explicitStartPos    = [];       % Explicit start position
centerTargetState   = [];       % State containing the target to use as the workspace center
binWidth            = 45;       % Bin width for binning spikes (ms)
spikeOnsetLag       = 500;      % Pre-onset time included in binned spike data (ms)
primaryCursor       = 'Neural Decoder';
verbose             = true;     % Display messages

% Handle case where only the Trial data is provided as an input.  In this
% case, use the 'getTaskProcessingParams' function to get the required
% inputs.
if nargin == 1
    fprintf('Automatically determining processing parameters from trial name.\n')
    P = getTaskProcessingParams(T(1).Overview.trialName);
    stateNames = P.tdParams.stateNames;
    targetState = P.tdParams.targetState;
    successOnly = P.tdParams.successOnly;
    rmNoSpikeData = P.tdParams.rmNoSpikeData;
    
    % Get center target state if specified
    if isfield(P.tdParams,'centerTargetState')
        centerTargetState = P.tdParams.centerTargetState;
    end
    varargin = P.tdParams.varargStr;
end

% Parse optional agruments
assignopts(who,varargin);

fprintf('Converting to TrajectoryData object ... ')

% If an explicit start position is defined, use this to re-center the
% workspace.
useExplicitStartPos = false;
if ~isempty(explicitStartPos)
    useExplicitStartPos = true;
    disp('Using explicit start position ...')
end

% If the center target state is defined, use the target for this state to
% re-center the position data in the workspace.
centerStateDefined = false;
if ~isempty(centerTargetState)
    centerStateDefined = true;
    fprintf('Using start position during state ''%s'' ...',centerTargetState)
end

% If an explicit start position is not provided AND the center target state
% has not been specified, then calculate the workspace center as the
% average over all unique targets.
if ~useExplicitStartPos && ~centerStateDefined
    fprintf('Calculating start position from target positions ...')
end

% Remove invalid trials
T = removeInvalidTrials(T,successOnly,rmNoSpikeData);
nTrials = length(T);

% Adjust pre-onset bin time so that t = 0 will be a bin edge
spikeOnsetLag = floor(spikeOnsetLag/binWidth) * binWidth;
% For offline decoding purposes, there should be a lag used due to neural
% activity leading kinematics by ~100ms.  In order to allow this to be
% accounted for, the 'spikeOnsetLag' is used to include spike counts before
% the onset of the specified state.

% Parse state names.  targStateName can be a cell array to get data from 
% multiple states.In this case, the first entry is the start state and the
% last state is the end state.
if ~iscell(stateNames) % Convert to a cell array if a string is provided
    stateNames = {stateNames};
end
nStates = length(stateNames);

% Verify that the desired states are present for all trials.
onsetFlag = false(nTrials,1);
offsetFlag = false(nTrials,1);
trialStOnset = cell(nTrials,1);
trialStOffset = cell(nTrials,1);
statesFound = cell(nTrials,1);
for i = 1:nTrials
    % Need to define a target state index -- this is the state containing
    % the target for the trajectory.  It is important to ensure this is
    % present for *all* trials, not just successful ones.
    
    % Loop over all states names and get onset and offset.  The
    % 'getStateTiming' function has been updated to get the last instance
    % of the states of interest.  This is because we sometimes use looping
    % states to allow the animal to re-acquire the target during the hold
    % period.  In this case, we want to grab the last occurance of the
    % state.
    allStOnset = nan(nStates,1);
    allStOffset = nan(nStates,1);
    for j = 1:nStates
        [tempOnset,tempOffset] = getStateTiming(T(i),stateNames{j}, ...
            'mode','last');
        if ~isempty(tempOnset)
            allStOnset(j) = tempOnset;
        end
        if ~isempty(tempOffset)
        	allStOffset(j) = tempOffset;
        end
    end
    statesFound{i} = (~isnan(allStOnset)) & (~isnan(allStOffset));
    allStOnset = allStOnset(statesFound{i});
    allStOffset = allStOffset(statesFound{i});
    
    % Put state onset and offset into cell array.  If at least one of the
    % desired states is found, the trial will be included.  Otherwise, it
    % will be excluded from the TrajectoryData object.
    if ~isempty(allStOnset)
        trialStOnset{i} = allStOnset;
        onsetFlag(i) = true;
    end
    
    if ~isempty(allStOffset)
        trialStOffset{i} = allStOffset;
        offsetFlag(i) = true;
    end
end

% Need to cut out trials where the desired onset and offset weren't found
stateMask = (onsetFlag & offsetFlag);
if sum(stateMask) ~= length(stateMask)
    nInvalidTrials = sum(~stateMask);
    warning('%d invalid trial(s) found. These will be excluded.', ...
        nInvalidTrials);
    
    % Get rid of invalid trials
    T = T(stateMask);
    trialStOnset = trialStOnset(stateMask);
    trialStOffset = trialStOffset(stateMask);
    statesFound = statesFound(stateMask);
end

% Loop over trials and get trajectory data
nTrials = length(T);

% Create empty trajectory data object array
TD = repmat(TrajectoryData,nTrials,1);

for i = 1:nTrials
    % Get state onset and offset data from cell array
    allStOnset = trialStOnset{i};
    allStOffset = trialStOffset{i};
    
    % Currently, there is no check to verify that all intermediate states
    % are found.  The onset and offset for the trajectory are defined as
    % the onset for the first state found and the offset for the last state
    % found.  In the future this could be updated such that the states
    % found are concatenated together, even if they are not continuous.
    
    % Define onset and offset
    onset = allStOnset(1);
    offset = allStOffset(end);
    
    % Get target state index and name
    endStIdx = findStateIndex(T(i).Parameters.states,targetState);
    endTargetName = T(i).Parameters.states(endStIdx).Cursor.name;

    % Get movement onset and offset
    moveOnset = T(i).Events.timeMoveOnset;
    moveOffset = T(i).Events.timeMoveEnd;
    
    % Determine cursors.  Some tasks will have both a phasespace and a hand
    % cursor.  If this is the case, choose the appropriate cursor based on
    % the 'primaryCursor' optional argument.  Currently this defaults to
    % the Neural Decoder.
    cursors = {T(i).Parameters.states(endStIdx).Cursor.marker};
    nCursor = length(cursors);
    cursorIdx = 1;
    if nCursor > 1
        cursorIdx = find(strcmp(cursors,primaryCursor));
        % In case the primaryCursor index isn't one of the cursors reset
        % cursorIdx to 1
        if isempty(cursorIdx)
            fprintf(['Unable to find the primary Cursor: ' primaryCursor ...
                '... Setting the cursor index to 1'])
            cursorIdx = 1;
        end
    end
    
    % Get target data
    [tPos,tSz] = getTargetPosition(T(i),targetState,endTargetName, ...
        'cursorIndex',cursorIdx);  
    
    % Get start position from center target state if specified
    startPosDefined = false;
    sPos = [];
    if centerStateDefined
        % Get target name
        ctrStIdx = findStateIndex(T(i).Parameters.states,centerTargetState);
        centerTargetName = T(i).Parameters.states(ctrStIdx).Cursor(cursorIdx).name;
        % Get target position
        [sPos,~] = getTargetPosition(T(i),centerTargetState,centerTargetName, ...
            'cursorIndex',cursorIdx);
        startPosDefined = true;
    end

    if useExplicitStartPos
        sPos = explicitStartPos;
        % Invert explicit y-position
        if T(i).Overview.invertYAxis
            sPos(:,2) = -sPos(:,2);
        end
        startPosDefined = true;
    end
    
    % Get intermediate targets
    tempStatesFound = statesFound{i};
    tempStates = stateNames(tempStatesFound);
    intTargPos = nan(length(tPos),sum(tempStatesFound));
    intTargSz = nan(1,sum(tempStatesFound));
    for j = 1:length(tempStates)
        % Get intermediate target state and name
        intStateName = tempStates{j};
        intStateIdx = findStateIndex(T(i).Parameters.states,intStateName);
        intTargetName = T(i).Parameters.states(intStateIdx).Cursor(cursorIdx).name;
        
        % Get intermediate target pos
        [intTargPos(:,j),intTargSz(j)] = getTargetPosition(T(i), ...
            intStateName,intTargetName,'cursorIndex',cursorIdx);
    end
        
    % Get kinematic data.  Loop through all possible kinematic sources
    % (Phasespace, Neural Decoder, and Auto-monkey) and get kinematic data
    % if available.  This allows both hand and brain control kinematics to
    % be plotted.
    
    % Determine the source of the kinematic data used to progress through
    % the state.
    kinSource = T(i).Parameters.states(endStIdx).Cursor(cursorIdx).marker;
    kinData{1} = T(i).Data.kinematics.cursor;
    kinDataSource{1} = 'hand';
    switch kinSource
        case 'Phasespace'
            % Don't need to do anything else
        case 'Neural Decoder'
            % Get decoder output in addition to hand position
            kinData{2} = T(i).Data.decoder.cursor;
            kinDataSource{2} = 'brain';
        case 'Auto-monkey'
            % Get auto-monkey position in addition to hand position
            kinData{2} = T(i).Data.autoMonkey;
            kinDataSource{2} = 'brain';
        case 'Force Cursor'
            % Get data from force object
            kinData{2} = T(i).Data.kinematics.force;
            kinDataSource{2} = 'force';
    end
    
    % Loop over kinematic sources and get data.  This could probably be
    % simplified by creating a vector of 'KinematicData' objects.
    %
    % Note -- currently this does not get the *actual* force information
    % for the ForceData object -- it only gets the force cursor position.
    % Updating the function to do this will involve re-thinking how this is
    % handled, as the force cursor and force data have different sampling
    % rates (the force cursor is updated every ~3ms, while the force data
    % is sampled at 1kHz).
    for j = 1:length(kinData)    
        % Create time mask
        t = kinData{j}.time;
        tMask = (t >= onset) & (t <= offset);

        % Hack to get one additional sample - it appears that there is
        % something odd with the state transition times
        timingFix = 0;
        if timingFix
            dT = diff(tMask);
            offsetInd = find(dT < 0,1,'first');
            offsetInd = offsetInd + 1;
            tMask(offsetInd) = 1;
        end

        % Get trajectory data - currently only position and velocity
        p = kinData{j}.pos(tMask,:);
        v = kinData{j}.vel(tMask(1:size(kinData{j}.vel,1)),:); % Note: this could throw an error (could be different dimensionality than position)
        a = kinData{j}.acc(tMask(1:size(kinData{j}.acc,1)),:);
        t = t(tMask);
    
        % Put data into KinematicData object
        srcStr = [kinDataSource{j} 'Kin'];
        TD(i).(srcStr).time = t;
        TD(i).(srcStr).pos = p;
        TD(i).(srcStr).vel = v;
        TD(i).(srcStr).acc = a;
        TD(i).(srcStr).source = kinDataSource{j};
    end
    
    % Get neural data and bin.  The bin times used here will be in the same
    % time units as the onset/offset times.
    S = T(i).Data.spikes;  % Structure array of spike times
    w = [(onset-spikeOnsetLag) offset];
    [S,spikeCounts,spikeTime] = S.binSpikes(binWidth,onset,0,w);
    
    % Put data into trajectory data object
    TD(i).states = tempStates;
    TD(i).stateOnset = allStOnset;
    TD(i).stateOffset = allStOffset;
    TD(i).controlSource = kinSource;
    TD(i).kinTime = t;          % Currently the last source
    TD(i).pos = p;              % Currently the last source
    TD(i).vel = v;              % Currently the last source
    TD(i).acc = a;              % Currently the last source
    TD(i).startPosDefined = startPosDefined;
    TD(i).startPos = sPos';
    TD(i).targPos = tPos';
    TD(i).targSize = tSz;
    TD(i).intTargPos = intTargPos;
    TD(i).intTargSz = intTargSz;
    TD(i).tube = T(i).Parameters.tube;
    %TD(i).cursorSize = cSz; % This shouldn't be necessary
    TD(i).spikes = S;
    TD(i).binTime = spikeTime;
    TD(i).binnedSpikes = spikeCounts;
    TD(i).binnedChannel = [S.channel];
    TD(i).binnedSort = [S.sort];
    
    % Get decode information.  This is somewhat redundant.
    if strcmp(kinSource,'Neural Decoder')
        TD(i).decodeTime = T(i).Data.decoder.rtTime(tMask);
        TD(i).decodeSpikeCounts = T(i).Data.decoder.rawSpikeBins(tMask,:);
        TD(i).decodeState = T(i).Data.decoder.rawDecode(tMask,5:end);
        TD(i).decoderName = T(i).Data.decoder.parameters.name;
        TD(i).decoderBinWidth = T(i).Data.decoder.parameters.binwidth;
        
        % Find decoder number.  Generally, this is appended to the end of
        % the decoder name string.
        decNameStr = TD(i).decoderName;
        if ~isempty(decNameStr) % Handle case where decoder name is empty
            decNumIdx = strfind(decNameStr,'_');
            decNumStr = decNameStr(decNumIdx(end)+1:end);
            TD(i).decoderNum = str2num(decNumStr);
        end
    end
    
    % Get events.
    TD(i).trajOnset = onset;
    TD(i).trajOffset = offset;
    TD(i).goCue = T(i).Events.timeGoCue;
    TD(i).moveOnset = moveOnset;
    TD(i).moveOffset = moveOffset;
    TD(i).targetAcquired = T(i).Events.timeTargetAcquired;
    TD(i).trialLen = T(i).Data.stateTransitions(2,end);
    
    % Get success and tag from trial data
    TD(i).trialID = str2num(T(i).Overview.trialNumber(6:end));
    TD(i).tag = T(i).Overview.tag;
    TD(i).successful = T(i).Overview.status;
    
    % Get meta data
    TD(i).subject = T(i).Overview.subject;
    TD(i).date = T(i).Overview.date;
    TD(i).trialName = T(i).Overview.trialName;
end

% Calculate movement onset for phasespace data - only do this if kinematic
% data is present
% KD = [TD.handKin];
% if ~isempty(KD)
%     KD = KD.calcMoveOnset;
%     TD = TD.setKinematicData(KD);
% end

% Print message that conversion is complete
fprintf(' done.\n')