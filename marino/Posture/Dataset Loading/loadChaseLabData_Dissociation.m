function [Data,postureIDs] = loadChaseLabData_Dissociation(dataset)

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = true;
        exclCh = [];
        exclZero = false; %Exclude channels that only contain sort 0
    getMarker = true;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    animalInitial = dataset(1);
    if strcmp(animalInitial,'N')
        trialName = 'Nigel Dissociation';
    else strcmp(animalInitial,'R')
       trialName = 'Rocky Dissociation'; 
    end

 %% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = '';
    postureIDs(2).ID = 2; postureIDs(2).Posture = '';
    
%% Load, Preprocess, and Label Data    
    date = dataset(2:end);
    month = date(5:6);
    if strcmp(animalInitial,'N')
        datasetDir = fullfile('D:\Animals\Nigel\2019',month,date);
    elseif strcmp(animalInitial,'R')
        datasetDir = fullfile('D:\Animals\Rocky\2020',month,date);
    end
    
    dataFileInfo = dir(fullfile(datasetDir,'*preprocessedData*'));
    dataPath = fullfile(datasetDir,dataFileInfo.name);
    rawData = load(dataPath);
    Data = getDataStructChaseLabDissociation20220419(rawData,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exclTrials',exclTrials);  

    
end