function [Data] = loadEarlData20210706_20211210

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
          %     [44 87 88 77 78 71 67 69 118]
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        checkPhasespaceSync = true;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = true;
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
    removeScreenFreezeTrials = false; %No photodiode data here, so you can't use it to determine whether screen froze
    inclStateTable = false;
    exclTrials = [];
 

%% Load, Preprocess, and Label Data 
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
   
    Data = getDataStruct20211210(Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
        'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
        'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
        'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);   
    
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    
    [Data,postureIDs] = labelPostures20210706(Data);
    
%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);

end