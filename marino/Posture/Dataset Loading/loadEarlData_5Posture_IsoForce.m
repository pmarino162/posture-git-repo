function [Data,postureIDs] = loadEarlData_5Posture_IsoForce(dataset)

%Loads and preprocesses Earl's 5 posture IsoForce datasets

%% Get that session exclude channels 
    switch dataset
        case 'E20200116'
            exclCh = [3 12 16 20 31 62 73 94];
        case 'E20200117'
            exclCh = [3 12 16 20 31 62 73 94];
        case 'E20200120'
            exclCh = [3 12 16 20 31 62 73 94];
    end
    
%% Get exclude trials for that session
    switch dataset
        case 'E20200116'
            exclTrials = [];
        case 'E20200117'
            exclTrials = [706];
        case 'E20200120'
            exclTrials = [];
    end

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = false;
        centerMarker = false; %Subracts workspace center from marker pose if true
    getForce = true;
        forceSetup = 'EL';
    getAlg = false;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = false;
    checkDecode = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;

%% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'I30';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'I15';
    postureIDs(3).ID = 3; postureIDs(3).Posture = 'N00';
    postureIDs(4).ID = 4; postureIDs(4).Posture = 'E15';
    postureIDs(5).ID = 5; postureIDs(5).Posture = 'E30';
    postureNames = {postureIDs.Posture};
    

%% Load and process each block of data (located in folders ending with 'grid'). Add to Data struct. Save GPFA and decoder params.
    Data = [];
    date = dataset(2:end);
    datasetDir = fullfile('D:\Animals\Earl\2020\01\',date);
    subFolderNames = dir(datasetDir); 
    block = 1; %Used for supplemental figures
    for i=6:10
        %Load data
        dataFileInfo = dir(fullfile(datasetDir,subFolderNames(i).name,'*translated*'));
        dataPath = fullfile(datasetDir,subFolderNames(i).name,dataFileInfo.name);
        tempData = load(dataPath);
        %Remove artifact
        [artifactTrials] = detectArtifact(tempData);
        tempData.Data(artifactTrials) = [];  
                
        %Get Data Struct
        tempData = getDataStruct20211210(tempData.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);      
        %Add Posture Info
        dataFileName = dataFileInfo.name;
        
        firstInd = strfind(dataFileName,'_isoForce') + 9;
        curInd = firstInd;
        while ~strcmp(dataFileName(curInd),'_')
            curInd = curInd + 1;
        end
        lastInd = curInd;          
        postureStr = dataFileName(firstInd:lastInd-1);
        
        if strcmpi(postureStr,'N') || strcmpi(postureStr,'neutral')
            postureStr = 'N00';
        end
        
        postureID = cellfun(@(x) strcmpi(x,postureStr),postureNames);
        postureID = postureIDs(postureID).ID;
        numTrials = size(tempData,2);
        for trial = 1:numTrials
            tempData(trial).conditionData.postureID = postureID;
            tempData(trial).conditionData.posture = postureStr;
            tempData(trial).conditionData.block = block;
        end

        %Add block of data to Data struct
        Data = [Data,tempData];
        block = block + 1;
    end
    
%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    
end