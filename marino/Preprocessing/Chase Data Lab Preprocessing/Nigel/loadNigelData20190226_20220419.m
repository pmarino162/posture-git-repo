function [Data,postureIDs] = loadNigelData20190226_20220419

 %% Preprocessing Parameters 
    dataset = 'N20190226';
    getSpikes = true;
        getSorts = true;
        exclCh = [];
        exclZero = false; %Exclude channels that only contain sort 0
    getMarker = true;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = 'Nigel Dissociation';


 %% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = '';
    postureIDs(2).ID = 2; postureIDs(2).Posture = '';
    
%% Load, Preprocess, and Label Data    
    rawData = load('D:\Animals\Nigel\2019\02\20190226\Nigel_20190226_preprocessedData_20220506_080540.mat');
    Data = getDataStructNigelReaching20220419(rawData,'trialName',trialName,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  

    
end