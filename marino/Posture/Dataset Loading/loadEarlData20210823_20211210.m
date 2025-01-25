function [Data] = loadEarlData20210823_20211210

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh = [7 42 43 47 51 58 65:78 81:87 89:91];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        %BC
        checkPhasespaceSyncBC = false;
        droppedMarkerTrialsBC = 'all';
        removeBadPhasespaceBC = false;
       centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = false;
        forceSetup = '';
    getAlg = false;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = true;
    %BC
    checkDecodeBC = true;
    %HC
    checkDecodeHC = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];
    
    
    %% Load, Preprocess, and Label Data 
    block1 = load('D:\Animals\Earl\2021\08\20210823\N00_brainControl\Earl20210823_N00_brainControl_SI_translated.mat');
    block1 = getDataStruct20211210(block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);        
    for trial = 1:size(block1,2)
        block1(trial).conditionData.posture = 'N00';
        block1(trial).conditionData.postureID = 1;
        block1(trial).conditionData.decoderID = 1;
    end
    
    block2 = load('D:\Animals\Earl\2021\08\20210823\N00_brainControl_02\Earl20210823_N00_brainControl_02_SI_translated.mat');
    block2 = getDataStruct20211210(block2.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);        
    for trial = 1:size(block2,2)
        block2(trial).conditionData.posture = 'N00';
        block2(trial).conditionData.postureID = 1;
        block2(trial).conditionData.decoderID = 2;
    end
    
    block3 =load('D:\Animals\Earl\2021\08\20210823\N00_brainControl_03\Earl20210823_N00_brainControl_03_SI_translated.mat');
    block3 = getDataStruct20211210(block3.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);        
    for trial = 1:size(block3,2)
        block3(trial).conditionData.posture = 'N00';
        block3(trial).conditionData.postureID = 1;
        block3(trial).conditionData.decoderID = 3;
    end
    
    block4 = load('D:\Animals\Earl\2021\08\20210823\N00_brainControl_04\Earl20210823_N00_brainControl_04_SI_translated.mat');
    block4 = getDataStruct20211210(block4.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);        
    for trial = 1:size(block4,2)
        block4(trial).conditionData.posture = 'N00';
        block4(trial).conditionData.postureID = 1;
        block4(trial).conditionData.decoderID = 4;
    end
    
    %% Combine data; keep only successful
    Data = [block1 block2 block3 block4];
    
end