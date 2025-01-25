function [Data,GPFAParams,decoderParams,postureIDs] = loadEarlData_20190830

%Loads and preprocesses Earl's 5 posture BCI CO datasets

%% Get that session exclude channels 
    dataset = 'E20190830';
    switch dataset
        case 'E20190830'
            exclCh = [];
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
    getKin = true;
    checkDecode = true;
    removeScreenFreezeTrials = true;
    inclStateTable = false;
    exclTrials = [];

%% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'I50';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'I20';
    postureIDs(3).ID = 3; postureIDs(3).Posture = 'N00';
    postureIDs(4).ID = 4; postureIDs(4).Posture = 'E20';
    postureIDs(5).ID = 5; postureIDs(5).Posture = 'E50';
    postureNames = {postureIDs.Posture};
    
%% Set up decoderParams and GPFAParams
    decoderParams = struct('postureID',[],'params',[]);
    GPFAParams = struct('postureID',[],'params',[]);
    paramsInd = 1;
    
%% Load and process each block of data (located in folders ending with 'grid'). Add to Data struct. Save GPFA and decoder params.
    Data = [];
    date = dataset(2:end);
    datasetDir = fullfile('D:\Animals\Earl\2019\08\',date);
    subFolderNames = dir(datasetDir); 
    block = 1; %Used for supplemental figures
    for i=1:numel(subFolderNames)
        %Get data from folders ending with 'CO'
        endsWithCO = endsWith(subFolderNames(i).name,'CO');
        if endsWithCO
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
            lastInd = strfind(dataFileName,'CO')-1;
            curInd = lastInd-1;
            while ~strcmp(dataFileName(curInd),'_')
                curInd = curInd -1;
            end
            firstInd = curInd + 1;
            postureStr = dataFileName(firstInd:lastInd);
            if strcmpi(postureStr,'neutral')
                postureStr = 'N00';
            elseif strcmpi(postureStr,'internal20')
                postureStr = 'I20';
            elseif strcmpi(postureStr,'internal50')
                postureStr = 'I50';
            elseif strcmpi(postureStr,'external20')
                postureStr = 'E20';
            elseif strcmpi(postureStr,'external50')
                postureStr = 'E50';
            end
            postureID = cellfun(@(x) strcmpi(x,postureStr),postureNames);
            postureID = postureIDs(postureID).ID;
            numTrials = size(tempData,2);
            for trial = 1:numTrials
                tempData(trial).conditionData.postureID = postureID;
                tempData(trial).conditionData.posture = postureStr;
                tempData(trial).conditionData.block = block;
            end
            %Load GPFA and Decoder Params; Save
            %I30GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run005\gpfa_xDim10.mat');
            blockNum = dataFileName(firstInd-3:firstInd-2);
            
            tempGPFAParams = load(fullfile([datasetDir,'\analysis\mat_results\run0',blockNum],'gpfa_xDim10.mat'));
            GPFAParams(paramsInd).postureID = postureID;
            GPFAParams(paramsInd).params = tempGPFAParams;
            
            tempDecoderParams = tempData(1).Decoder.Parameters;
            decoderParams(paramsInd).postureID = postureID;
            decoderParams(paramsInd).params = tempDecoderParams;
            paramsInd = paramsInd + 1;
            
            %Add block of data to Data struct
            Data = [Data,tempData];
            block = block + 1;
        end

    end
    
%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    
end