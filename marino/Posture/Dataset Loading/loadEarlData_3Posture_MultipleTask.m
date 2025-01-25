function [Data] = loadEarlData_3Posture_MultipleTask(dataset)

%Loads and preprocesses Earl's 3-posture Multiple Task datasets

%% Get exclude channels for each session
    switch dataset
        case 'E20200311'
            exclCh = [3 12 16 20 31 73 92 94];
        case 'E20200312'
            exclCh = [3 12 16 20 31 73 92 94];
        case 'E20200313'
            exclCh = [3 12 16 20 31 73 92 94];
        case 'E20200314'
            exclCh = [3 12 16 20 31 73 85 92 94];
    end

%% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        %BC
        checkPhasespaceSyncBC = false;
        droppedMarkerTrialsBC = 'all';
        removeBadPhasespaceBC = false;
        %HC
        checkPhasespaceSyncHC = true;
        droppedMarkerTrialsHC = 'no drop';
        removeBadPhasespaceHC = true;
        %Iso
        checkPhasespaceSyncIso = false;
        droppedMarkerTrialsIso = 'all';
        removeBadPhasespaceIso = false;
        centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = true;
        forceSetup = 'EL';
    getAlg = true;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = true;
    %BC
    checkDecodeBC = true;
    %HC
    checkDecodeHC = false;
    %Iso
    checkDecodeIso = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];

%% Load and preprocess each block of data
    Data = [];
    date = dataset(2:end);
    datasetDir = fullfile('D:\Animals\Earl\2020\03\',date);
    subFolderNames = dir(datasetDir); 
    block = 1; %Used for supplemental figures
    for i=1:numel(subFolderNames)
        %Get block, posture, and task
        undInd = strfind(subFolderNames(i).name,'_');
        task = '';
        if length(undInd) > 1
            posture = subFolderNames(i).name(undInd(1)+1:undInd(2)-1);
            task = subFolderNames(i).name(undInd(2)+1:end);
        end
        %If the task is one that we're analyzing, process this block
            %(Except for neutral forceBar_CO on 2020311, bc I later re-ran
            % this block with updated target locations)
        if ismember(task,{'grid','forceBar_CO','isoForce','forceBar_CO_02'}) && ~(strcmpi(date,'20200311') && strcmpi(posture,'neutral') && strcmpi(task,'forceBar_CO'))
            %Load data 
            if strcmpi(date,'20200311')
                dataFileInfo = dir(fullfile(datasetDir,subFolderNames(i).name,'*SI_UDP_translated*'));
            else
                dataFileInfo = dir(fullfile(datasetDir,subFolderNames(i).name,'*translated*'));
            end
            dataPath = fullfile(datasetDir,subFolderNames(i).name,dataFileInfo.name);
            tempData = load(dataPath);
            %Remove artifact
            [artifactTrials] = detectArtifact(tempData);
            tempData.Data(artifactTrials) = []; 
            %Remove problem trials (before preprocessing)
            switch dataset
                case 'E20200311'
                    [tempData] = cleanData20200311(tempData);
                case 'E20200312'
                    [tempData] = cleanData20200312(tempData);
                case 'E20200313'
                case 'E20200314'
                    [tempData] = cleanData20200314(tempData);
            end
            
            %Get Data Struct
            switch task
                case 'grid'
                    checkPhasespaceSync = checkPhasespaceSyncBC;
                    droppedMarkerTrials = droppedMarkerTrialsBC;
                    removeBadPhasespace = removeBadPhasespaceBC;
                    checkDecode = checkDecodeBC;
                case {'forceBar_CO','forceBar_CO_02'}
                    %Rename trials
                    numTrials = size(tempData.Data,2);
                    for trial = 1:numTrials
                        tempData.Data(trial).Overview.trialName = 'HC_CenterOut_ForceBar_20200314';
                    end
                    checkPhasespaceSync = checkPhasespaceSyncHC;
                    droppedMarkerTrials = droppedMarkerTrialsHC;
                    removeBadPhasespace = removeBadPhasespaceHC;
                    checkDecode = checkDecodeHC;
                case 'isoForce'
                    checkPhasespaceSync = checkPhasespaceSyncIso;
                    droppedMarkerTrials = droppedMarkerTrialsIso;
                    removeBadPhasespace = removeBadPhasespaceIso;
                    checkDecode = checkDecodeIso;
            end
            tempData = getDataStruct20211210(tempData.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
                'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
                'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
            %Add posture and task info
            switch posture
                case 'neutral'
                    postureStr = 'N00';
                    postureID = 2;
                case 'external'
                    postureStr = 'E30';
                    postureID = 3;
                case 'internal'
                    postureStr = 'I30';
                    postureID = 1;
            end
            switch task
                case 'grid'
                    taskStr = 'BC';
                    taskID = 1;        
                case {'forceBar_CO','forceBar_CO_02'}
                    taskStr = 'HC';
                    taskID = 2;
                case 'isoForce'
                    taskStr = 'Iso';
                    taskID = 3;
            end            
            numTrials = size(tempData,2);
            for trial = 1:numTrials
                tempData(trial).conditionData.postureID = postureID;
                tempData(trial).conditionData.posture = postureStr;
                tempData(trial).conditionData.task = taskStr;
                tempData(trial).conditionData.taskID = taskID;
                tempData(trial).conditionData.block = block;
                %For BCI task, rename target1ID targetID
                if strcmpi(task,'grid')
                    tempData(trial).targetData.targetID = tempData(trial).targetData.target1ID;
                end
            end

            %Add block of data to Data struct
            Data = [Data,tempData];
            block = block + 1;
        end
    end

%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    
end