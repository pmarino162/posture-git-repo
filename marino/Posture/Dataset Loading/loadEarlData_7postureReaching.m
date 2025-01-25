function [Data] = loadEarlData_7postureReaching(dataset)

%% Get that session exclude channels 
    switch dataset
        case 'E20210706'
            exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
        case 'E20210707'
            exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
        case 'E20210708'
            exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
        case 'E20210709'
            exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
        case 'E20210710'
            exclCh = [25 26 44 50 53 67 69 73 78 79 87 88 100 116 117 118 120 123];
    end
    
 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
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
 
%% Load dataset
    switch dataset
        case 'E20210706'
            load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
        case 'E20210707'
            load('D:\Animals\Earl\2021\07\20210707\Earl20210707_gridReaching_SI_translated.mat');
        case 'E20210708'
            load('D:\Animals\Earl\2021\07\20210708\Earl20210708_gridReaching_SI_translated.mat');
        case 'E20210709'
            load('D:\Animals\Earl\2021\07\20210709\Earl20210709_gridReaching_SI_translated.mat');
        case 'E20210710'
            load('D:\Animals\Earl\2021\07\20210710\Earl20210710_gridReaching_SI_translated.mat');
    end
    
%% Preprocess Data  
    Data = getDataStruct20211210(Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
        'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
        'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
        'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);   

%% Clean data
    Data = cleanEarl_7postureReachingData(Data,dataset);
    [Data] = cleanData20210706(Data);

%% Label data
    [Data,postureIDs] = labelPostures20210706(Data);
    
%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);

end