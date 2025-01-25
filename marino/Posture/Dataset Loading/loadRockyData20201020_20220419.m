function [Data,postureIDs] = loadRockyData20201020_20220419

 %% Preprocessing Parameters 
    dataset = 'R20201020';
    getSpikes = true;
        getSorts = true;
        exclCh = [];
        exclZero = false; %Exclude channels that only contain sort 0
    getMarker = false;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = 'Rocky Posture BC Center Out';

 %% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'N00';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'A90';
    
%% Load, Preprocess, and Label Data    
    Data = [];
    %Posture 1 (block 1)
        %IntData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390_intermediate.mat');
        EZData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390.raw_EZ.mat');
        tempData = getDataStructRockyBC_EZ_20220419(EZData.EZ,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exclTrials',exclTrials);  
         
%         rawData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390.mat');
%         tempData = getDataStructRockyBC20220419(rawData.Data,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
%             'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
        
        %Remove first 40 trials (calibration)
        tempData(1:40) = [];
        %Add Posture
        for i=1:size(tempData,2)
            tempData(i).conditionData.postureID = 1;
            tempData(i).conditionData.posture = 'N00';
        end
        %Add to Data Struct
        Data = [Data,tempData];
        

    %Posture 2 (block 2)
%         rawData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00393.mat');
%         tempData = getDataStructRockyBC20220419(rawData.Data,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
%             'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
        EZData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00393.raw_EZ.mat');
        tempData = getDataStructRockyBC_EZ_20220419(EZData.EZ,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
        
        %Remove first 40 trials (calibration)
        tempData(1:40) = [];
        %Add Posture
        for i=1:size(tempData,2)
            tempData(i).conditionData.postureID = 2;
            tempData(i).conditionData.posture = 'A90';
        end
        %Add to Data Struct
        Data = [Data,tempData];
    
  %% Clean Data
    [Data] = cleanRockyBCIData(Data,dataset);
    
end