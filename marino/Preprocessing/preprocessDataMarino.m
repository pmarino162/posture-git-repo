function [Data,P] = preprocessDataMarino(Data, varargin)
% [Data] = preprocessData(Data, varargin)
%   Inputs: raw data structure, variable arguments
%   Outputs: preprocessed data structure
%
%   Variable Arguments:
%       SPIKES: true (default) --> preprocess neural activity
%               false --> do not preprocess neural activiity
%
%       GETUDP: true --> get spikes transmitted via GETUDP
%               false (default) --> get spikes saved on TDT computer
%
%       USESORTS: true (default) --> separate spiking data by channel and
%                                    by sort
%                 false --> separate data by channel only
%
%       SORTCODES: array of sort codes to get for analysis. default is
%                     [1:30].  0 corresponds to unsorted spikes, and 31
%                     corresponds to outliers.
%
%       PSPSAMPLEDIFFTHRESHOLD: this is the threshold (in # of samples) above which a
%                             time gap in phasespace data is considered to
%                             be a dropped sample. default = 1 sample.
%       DECODE: true (default) --> include the decoder information in the
%                   processed struct
%               false --> don't include the decoder information
%
%       DROPPEDMARKERTRIALS: all (default) --> default is to get all trials
%                            drop --> get only trials w/ dropped marker
%                            nodrop --> get only trials w/o dropped marker
%
%       SINGLEMARKER: true (default) --> analyze only the marker that is
%                       used for hand data
%                     false --> use all marker (code for this has not been
%                       developed yet)
%
%       DIMENSIONS: 3 (default) --> analyze three dimensions on hand data
%                   2 --> analyze two dimensions (x,y) of hand data
%
%       CHECKPHASESPACESYNC: true (default) --> remove trials for which
%                                                initial phasespace time is
%                                                out of sync
%                            false --> don't check phasespace sync

%       PHOTOTHRESHOLD: the voltage threshold used to ensure the diode(s)
%                           was working (default = 3.5V)
%
%       WAVEFORMS: true --> get the snippet waveforms
%                  false (default) --> don't get the snippet waveforms
%
%       DETERMINEATTEMPT: true (default) --> determine which trials were
%                                               attmpted by the subject
%                         false --> don't determine which trials were
%                                       attempted
%
%       DISPTEXT: true (default) --> display text while running preprocess
%                 false --> don't display text
%
%       REMOVESCREENFREEZETRIALS: true --> remove the trials on which the 
%                                              screen froze
%                                 false (default) --> don't remove the
%                                              screen freeze trials
%
%       ADDITIONALFREEZETRIALS: sometimes the trial folling a screen freeze
%                   trial is messed up. This is the number of trials
%                   following the screen freeze trials to remove. (default
%                   = 1).
%
%       AUTOMONKEY: true (default) --> default is to process the
%                       auto-monkey kinematics
%                   false --> don't process the auto-monkey kinematics
%
%       INTERPAUTOMONKEY: true --> interpolate the auto-monkey positions
%                         false (default) --> don't interpolate them
%
%       CHECKDECODE: true (default) --> check to make sure the timestamps
%                                       of the decoder are monotonically
%                                       increasing
%                    false --> don't do this check
%       EVENTPLOT: true --> creates plots of the difference between the
%                               phototransition and state event times
%                  false (default) --> doesn't create the plots.
%
% Updated from preprocessPTS by Patrick Sadtler.  This version determines
% the appropriate preprocessing parameters depending on the trial type.
% 
% Author:       Alan D. Degenhart
% Date Created: N/A
% Last Updates: 
% 2016.06.07
% 2017.03.05        EMG         Added PRIMARYMOVEONSETMODE option to computeMoveTimes
% 2019.04.26        NTM         Added ADAPTIVETHRESHOLD and
%                               ADAPTIVETHRESHOLDVALUE options to
%                               photoTransitions

% Set the default values for the optional arguments
SPIKES = true;                          % default is to preprocess neural activity
USESORTS = true;                        % default is to separate neural data by channel and sort
GETUDP = false;                         % default is to get spikes saved on TDT computer
TUBE = true;                            % default is to include the tube
SORTCODES = [-1:31];                    % default is to get spikes with a sort code in the interval [-1:31]
PSPSAMPLEDIFFTHRESHOLD = 30;             % default threshold for phasespace dropped trials
DROPPEDMARKERTRIALS = 'all';            % default is to get all trials
DECODE = true;                          % default is to include the decoder information in the processed struct
SINGLEMARKER = true;                    % default is to use only one marker when analyzing the hand data
DIMENSIONS = 3;                         % default is to analyze three dimensions of hand data
CHECKPHASESPACESYNC = true;             % default is to remove trials that are out of sync
PHOTOTHRESHOLD = 1.75;                  % default voltage to determine if photodiode(s) was working is 3.5V. Now 1.75V
ADAPTIVETHRESHOLD = true;               % default is to use adaptive threshold to detect photodiode state transitions
ADAPTIVETHRESHOLDVALUE = 0.5;           % default is to use 0.5 of max voltage as photodiode adaptive threshold value
WAVEFORMS = false;                      % default is to not get the snippet waveforms
DETERMINEATTEMPT = true;                % default is to determine which trials were attemptedby subject
DISPTEXT = true;                        % default is to display text when preprocessing
REMOVESCREENFREEZETRIALS = false;       % default is to remove the trials on which the screen froze
ADDITIONALFREEZETRIALS = 1;             % default is to remove the trial that occured after the screen froze
AUTOMONKEY = true;                      % default is to process the auto-monkey kinematics
INTERPAUTOMONKEY = true;                % default is to interpolate the auto-monkey positions
CHECKDECODE = true;                     % default is to check that the timestamps on the decoder are monotonically increasing
EVENTPLOT = false;                      % default is to not create plots.

varRemain = assignopts('ignorecase', who, varargin);

if ~isempty(varRemain)
    fprintf('\n\n')
    disp('The following optional arguments were invalid (preprocessData.m):')
    disp(varRemain)
    disp('Press any key to continue')
    pause
end

if DISPTEXT
    fprintf('\r')
end

% Check to see if waveform data exists.  If it does, include it.
nTrials = length(Data);
includesWaveforms = false(nTrials,1);
for i = 1:nTrials
    if ~isempty(Data(i).TrialData.TDT.snippetWaveforms)
        includesWaveforms(i) = true;
    end
end
if sum(includesWaveforms) > 0
    WAVEFORMS = true;
    fprintf('Including snippet waveforms.\n')
end

% Get preprocessing parameters.  These are trial-dependent parameters
% parameters which specify how the data is to be preprocessed.
nTrials = length(Data);
Ptemp = genPreprocessParams(Data(1));
P = repmat(Ptemp,nTrials,1);
for i = 1:nTrials
    P(i) = genPreprocessParams(Data(i));
end

% Check for trials where the phasespace data is out of sync.  Attempt to
% re-sync trials if possible.  Trials that cannot be re-synced will be
% discarded.
if CHECKPHASESPACESYNC
    fprintf('Checking data for phasespace sync ... ');
    Data = syncPSP(Data,P,'DISPTEXT',DISPTEXT);
    
    outOfSyncTrls = [];
    for syncTrl = 1:length(Data)
        if isfield(Data(syncTrl).Overview,'outOfSync')
            if Data(syncTrl).Overview.outOfSync && ~Data(syncTrl).Overview.resynced
                outOfSyncTrls = [outOfSyncTrls syncTrl];
            end
        end
    end
    Data = Data(setdiff([1:length(Data)],outOfSyncTrls));
    fprintf('done.\n')
end

% Check for trials where Phasespace marker data was dropped.
fprintf('Checking data for dropped markers ... ');
droppedMarkerTrials = findDroppedMarkers(Data,PSPSAMPLEDIFFTHRESHOLD, ...
    'SINGLEMARKER',SINGLEMARKER,'DISPTEXT',DISPTEXT);
for m = 1:length(Data)
    if ismember(m,droppedMarkerTrials)
        Data(m).Overview.droppedMarker = 1;
    else
        Data(m).Overview.droppedMarker = 0;
    end
end
fprintf('done.\n')

% Remove trials depending on the value of 'DROPPEDMARKERTRIALS'
if strcmp(DROPPEDMARKERTRIALS,'no drop')
    AllTrials=length(Data);
    Data = Data(setdiff([1:AllTrials],droppedMarkerTrials));
    P = P(setdiff([1:AllTrials],droppedMarkerTrials));
    disp([num2str(length(droppedMarkerTrials)),' trials with dropped markers removed from data.']);
elseif strcmp(DROPPEDMARKERTRIALS,'drop')
    Data = Data(droppedMarkerTrials);
    P = P(droppedMarkerTrials);
    disp([num2str(length(setdiff([1:length(Data)],droppedMarkerTrials))),' trials with no dropped markers removed data.']);
end

% % Determine if trials were attempted.
% if DETERMINEATTEMPT
%     fprintf('Finding attempted trials ... ');
%     Data = determineAttemptTrials(Data,P,'DISPTEXT',DISPTEXT);
%     fprintf('done.\n')
% end

% Verify that the online decode was correctly translated.  Discard trials
% where this is not the case.
if CHECKDECODE
    fprintf('Finding invalid decode trials ... ')
    invalidDecodeTrials = checkDecodeValidity(Data,'DISPTEXT',0);
    if ~isempty(invalidDecodeTrials)
        fprintf('%d trials with invalid decodes removed ... ',length(invalidDecodeTrials))
        Data(invalidDecodeTrials) = [];
        P(invalidDecodeTrials) = [];
    end
    fprintf('done.\n')
end
    
% Process neural data
if ~isempty(Data)
    % check to see if the neural data has already been processed
    if (~isfield([Data(1).TrialData],'spikes')) && SPIKES
        fprintf('Analyzing neural activity ... ');
        try
        Data = findUnits(Data,'USESORTS',true,'GETUDP',GETUDP,'SORTCODES',SORTCODES,...
            'WAVEFORMS',WAVEFORMS,'DISPTEXT',DISPTEXT);
        catch ERR
            keyboard
        end
        fprintf('done.\n')
    end
end

% Get LFP data
Data = getLFPData(Data);

% Photodiode data
if ~isempty(Data)
    if ~isfield([Data(1).TrialData],'photoTransitions')
        % check if the photodiode signals have already been processed
        fprintf('Checking photodiode signals ... ');
        [Data] = photoTransitions(Data,'THRESHOLD',PHOTOTHRESHOLD,...
            'DISPTEXT',DISPTEXT,'ADAPTIVETHRESHOLD',ADAPTIVETHRESHOLD,...
            'ADAPTIVETHRESHOLDVALUE',ADAPTIVETHRESHOLDVALUE);
        fprintf('done. \n');
    end
end

% Get event times
if ~isempty(Data)
    fprintf('Getting event times ... ');
    [Data,invalidTrials] = eventTimes(Data,P,'DISPTEXT',DISPTEXT,'PLOT',EVENTPLOT);
    fprintf('done.\n')
end

% Remove the trials during which the screen froze
if REMOVESCREENFREEZETRIALS
    % find the trials where the screen froze
    freezeTrials = findScreenFreezeTrials(Data,invalidTrials);
    if ~isempty(freezeTrials)
        % the trials immediately following a trial on which the screen
        % froze were usually bad (targets didn't show up in proper
        % locations), so those should be removed from the data
        extraFreezeTrials = freezeTrials + ADDITIONALFREEZETRIALS;
        allFreezeTrials = unique([freezeTrials extraFreezeTrials]);
        allFreezeTrials = find(allFreezeTrials <= length(Data));
        Data(allFreezeTrials) = [];
        P(allFreezeTrials) = [];
        fprintf('The screen froze on trials: ')
        fprintf('%d ',freezeTrials)
        fprintf('\rThose trials plus the folowing %d were moved from the data set.\r',ADDITIONALFREEZETRIALS)
    end
end

% Process decoded cursor kinematics
if ~isempty(Data)
    % check whether the decoded kinematics were already processed
    if (~isfield(Data(1).TrialData,'DecodeKinematics')) && DECODE
        fprintf('Getting decode kinematics ... ');
        Data = findBrainKinematics(Data,'VERBOSE',DISPTEXT);
        fprintf('done.\n')
    end
end

% Process auto monkey data
if ~isempty(Data)
    % check to make sure auto monkey data is present and has not been
    % processed
    if ~isfield(Data(1).TrialData,'AutoMonkeyKinematics') && AUTOMONKEY
        if isfield(Data(1).TrialData.Marker,'autoMonkey')
            if ~isempty(Data(1).TrialData.Marker.autoMonkey.position)
                fprintf('Getting auto-monkey kinematics ... ')
                Data = findAutoMonkeyKinematics(Data,'DIMENSIONS',DIMENSIONS, ...
                    'INTERP',INTERPAUTOMONKEY);
                fprintf('done.\n')
            end
        end
    end
end

% Add tube data for each state
if ~isempty(Data)
    Data = rearrangeTubeData(Data);
end

% Housekeeping for a few fields
for trial=1:length(Data)
    [Data(trial).Parameters.stateNames] = [Data(trial).Parameters.StateTable.stateName];
   
    Data(trial).Parameters.MarkerTargets = [Data(trial).Parameters.StateTable.Hand];

    % add field description
    Data(trial).Parameters.fieldDescription = {'Data.Parameters.TrialTargets indicates the targets there displayed on the screen to the subject.', ...
        'Data.Parameters.MarkerTargets indicates the targets to which the subject had to reach with the cursor.'};
end