function [trialsToUse] = subselectToRemoveArtifact(pN,fN,decoderID,calibFilename,trialsToUse)

% This function identifies a list of target-balanced, artifact-free, force bar failure-free trials
% for decoder calibration. 
%
% Inputs
%   - pN: path name for .mat file of translated data (with snippet info)
%   - fN: file name for .mat file of translated data (with snippet info)
%   - decoderID: the decoderID for the decoder that you'll calibrate with
%                the trials you keep
%   - calibFilename: the shortened file name for the translated .mat file used for
%                    calibration
%   - trialsToUse: a struct that stores trials to use for each decoder ID
%
% Outputs
%   - trialsToUse: a struct that stores trials to use for each decoder ID
%

%% Set threshold value and number of bins necessary for exlclusion
    threshold = 30;
    numBinsNec = 10;
    
%% Define state name for force bar failures
    forceBarFailureStateName = 'Failure Touch Bar';
    
%% Load data
    dataPath = fullfile(pN,fN);
    Data = load(dataPath);
    Data = Data.Data;

%% Preprocess data
    procData = preprocessData(Data);

%% Choose the set of trials to examine
    numTrials = size(procData,2);
    allTrials = [1:numTrials];
    trialsToExamine = allTrials;
    previousTrialsToUse = [];
    %If you've already calibrated with data this calibFile, only examine new trials. 
    if decoderID ~= 1
        previousCalibFilename = trialsToUse(decoderID-1).calibFilename;
        if strcmpi(calibFilename,previousCalibFilename)
            previousAllTrials = trialsToUse(decoderID-1).allTrials;
            trialsToExamine(previousAllTrials) = [];
            previousTrialsToUse = trialsToUse(decoderID-1).trialsToUse;
        end
    end
    numTrialsToExamine = size(trialsToExamine,2);   

%% Get number of channels and sorts; create channelSpikes struct
    spikes = procData(1).TrialData.spikes;
    channelList = unique([spikes.channel]);
    numChannels = length(channelList);
    numSorts = zeros(1,numChannels);
    for channel = 1:numChannels
        chSpikes = spikes([spikes.channel]==channel);
        numSorts(channel) = size(chSpikes,2);   
    end

%% Get target info for all trials
    for trial = allTrials
        trialName = procData(trial).Overview.trialName;
        TrialTargets = procData(trial).Parameters.TrialTargets;
        [targetData] = getTargetDataArtifactDetection(trialName,TrialTargets,Data);
        targetID = targetData.targetID;
        procData(trial).targetID = targetID;
    end
    
%% For each trial, check for artifact and forcebar failures.  Plot sum of spikes for each trial
    artifactTrials = [];
    forceBarFailures = [];
    f = figure;
    for trial = trialsToExamine
        %Get Binned Spikes
        spikes = procData(trial).TrialData.spikes;
        stateTransitions = procData(trial).TrialData.stateTransitions;
        [binTimes,allSpikeBins] = getBinnedSpikesArtifactDetection(spikes,stateTransitions,channelList,numChannels,numSorts);
        %Get targetID
        targetID = procData(trial).targetID;
        %Sum Across The Array 
        spikeSum = sum(double(allSpikeBins),2);
        %Plot Results
        subplot(ceil(sqrt(numTrialsToExamine)),ceil(sqrt(numTrialsToExamine)),trial-trialsToExamine(1)+1)
        plot(spikeSum)
        hold on
        ax = gca;
        plot(ax.XLim,[threshold,threshold],'r')
        xlabel('time (ms)')
        ylabel('# spikes')
        ylim([0 150])
        title(['Trial ',num2str(trial),newline,'Target ID: ',num2str(targetID)])
        %Compare to threshold value and flag trials containing artifact
        if sum(spikeSum > threshold) >= numBinsNec
            artifactTrials = [artifactTrials,trial];
        end
        %Flag trials containing force bar failures
        stateNames = procData(trial).Parameters.stateNames;
        stateTransitions = double(procData(trial).TrialData.stateTransitions(1,:));
        forceBarFailureStateInd = min(find([cellfun(@(x) strcmpi(x,forceBarFailureStateName),stateNames)]==1));
        if ismember(forceBarFailureStateInd,stateTransitions)
           forceBarFailures = [forceBarFailures,trial]; 
        end
    end
    f.Position = [0 0 1750 1000];
    
%% Create goodTrials, which have neither artifact nor force bar failures
    goodTrials = trialsToExamine;
    idx=ismember(goodTrials,[artifactTrials,forceBarFailures]);
    goodTrials(idx) = [];

%% Get Block Target Info
    targetIDs = [procData.targetID];
    possibleTargetIDs = [1:8];
    
%% Create bar plots of removals by target
    figure
    bar(possibleTargetIDs,histcounts(targetIDs(artifactTrials),[possibleTargetIDs,possibleTargetIDs(end)+1]))
    xlabel('Target ID')
    ylabel('Number of Removals')
    title('Artifact Removals By Target')
    
    figure
    bar(possibleTargetIDs,histcounts(targetIDs(forceBarFailures),[possibleTargetIDs,possibleTargetIDs(end)+1]))
    xlabel('Target ID')
    ylabel('Number of Removals')
    title('Force Bar Removals By Target')
    
%% Create targetBalancedTrials
    targetBalancedTrials = [];
    targetIDsGoodTrials = targetIDs(goodTrials);
    %Get minimum number of times any inidividual target appeared in
    %trials without artifact
    minNumber = numTrialsToExamine;
    for target = possibleTargetIDs
        numberOfOccurrences = double(sum(targetIDsGoodTrials==target));
        if numberOfOccurrences < minNumber
            minNumber = numberOfOccurrences;
        end
    end
    %If minimum number is > 4, set it equal to 4 (so we only get 32 trials
    %from the block)
    if minNumber > 4
        minNumber = 4;
    end
    %Choose First minNumber trials of each target ID
    for target = possibleTargetIDs
        validTrialInds = intersect(find(targetIDs==target),goodTrials);
        targetBalancedTrials = sort([targetBalancedTrials,validTrialInds(1:minNumber)]);
    end
    
%% Ask user whether to keep results; if so, populate trialsToUse
    %Print Resutls
    decoderID
    calibFilename
    artifactTrials
    forceBarFailures
    goodTrials
    numGoodTrials = size(goodTrials,2)
    targetBalancedTrials
    numTargetBalancedTrials = size(targetBalancedTrials,2)
    %Set flag to 1 when you're done talking with the user
    flag = 0;
    while flag == 0
        prompt = 'Save targetBalancedTrials for this decoderID? Y/N: ';
        str = input(prompt,'s');
        if strcmpi(str,'y')
            trialsToUse(decoderID).decoderID = decoderID;
            trialsToUse(decoderID).calibFilename = calibFilename;
            trialsToUse(decoderID).allTrials = allTrials;
            trialsToUse(decoderID).artifactTrials = artifactTrials;
            trialsToUse(decoderID).forceBarFailures = forceBarFailures;
            trialsToUse(decoderID).goodTrials = goodTrials;
            trialsToUse(decoderID).targetBalancedTrials = targetBalancedTrials;
            trialsToUse(decoderID).trialsToUse = [previousTrialsToUse,targetBalancedTrials];
            flag = 1;
        elseif strcmpi(str,'n')
            flag = 1;
        else  
        end
    end
    
end