function [Data,postureIDs] = loadEarlData20190729

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh =  [10 12 31];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = false;
        centerMarker = true; %Subracts workspace center from marker pose if true
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
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'E50';
    
%% Load, Preprocess, and Label Data   
    E50_E50d_CO_block1 = load('D:\Animals\Earl\2019\07\20190729\05_external_centerOut\Earl20190729_05_external_centerOut_SI_SW_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E50_E50d_CO_block1);
        E50_E50d_CO_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        E50_E50d_CO_block1 = getDataStruct20211210(E50_E50d_CO_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E50_E50d_CO_block1,2);
        for trial = 1:numTrials
           E50_E50d_CO_block1(trial).conditionData.postureID = 2;
           E50_E50d_CO_block1(trial).conditionData.posture = 'E50';
           E50_E50d_CO_block1(trial).conditionData.decoderPosture = 'E50';
           E50_E50d_CO_block1(trial).conditionData.decoderPostureID = 2;
        end

    I50_I50d_CO_block1 = load('D:\Animals\Earl\2019\07\20190729\10_internal_centerOut\Earl20190729_10_internal_centerOut_SI_SW_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I50_I50d_CO_block1);
        I50_I50d_CO_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        I50_I50d_CO_block1 = getDataStruct20211210(I50_I50d_CO_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I50_I50d_CO_block1,2);
        for trial = 1:numTrials
           I50_I50d_CO_block1(trial).conditionData.postureID = 1;
           I50_I50d_CO_block1(trial).conditionData.posture = 'I50';
           I50_I50d_CO_block1(trial).conditionData.decoderPosture = 'I50';
           I50_I50d_CO_block1(trial).conditionData.decoderPostureID = 1;
        end
    
    I50_I50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190729\12_internalPosture_internalDecoder_centerOutCenter_01\Earl20190729_12_internalPosture_internalDecoder_centerOutCenter_01_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I50_I50d_COC_block1);
        I50_I50d_COC_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        I50_I50d_COC_block1 = getDataStruct20211210(I50_I50d_COC_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I50_I50d_COC_block1,2);
        for trial = 1:numTrials
           I50_I50d_COC_block1(trial).conditionData.postureID = 1;
           I50_I50d_COC_block1(trial).conditionData.posture = 'I50';
           I50_I50d_COC_block1(trial).conditionData.decoderPosture = 'I50';
           I50_I50d_COC_block1(trial).conditionData.decoderPostureID = 1;
        end

    I50_E50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190729\13_internalPosture_externalDecoder_centerOutCenter\Earl20190729_13_internalPosture_externalDecoder_centerOutCenter_SI_SW_translated.mat');
    %Remove artifact
        [artifactTrials] = detectArtifact(I50_E50d_COC_block1);
        I50_E50d_COC_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        I50_E50d_COC_block1 = getDataStruct20211210(I50_E50d_COC_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I50_E50d_COC_block1,2);
        for trial = 1:numTrials
           I50_E50d_COC_block1(trial).conditionData.postureID = 1;
           I50_E50d_COC_block1(trial).conditionData.posture = 'I50';
           I50_E50d_COC_block1(trial).conditionData.decoderPosture = 'E50';
           I50_E50d_COC_block1(trial).conditionData.decoderPostureID = 2;
        end
        
    I50_I50d_COC_block2 = load('D:\Animals\Earl\2019\07\20190729\14_internalPosture_internalDecoder_centerOutCenter_02\Earl20190729_14_internalPosture_internalDecoder_centerOutCenter_02_SI_SW_translated.mat');
        [artifactTrials] = detectArtifact(I50_I50d_COC_block2);
        I50_I50d_COC_block2.Data(artifactTrials) = [];       
        %Get Data Struct
        I50_I50d_COC_block2 = getDataStruct20211210(I50_I50d_COC_block2.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I50_I50d_COC_block2,2);
        for trial = 1:numTrials
           I50_I50d_COC_block2(trial).conditionData.postureID = 1;
           I50_I50d_COC_block2(trial).conditionData.posture = 'I50';
           I50_I50d_COC_block2(trial).conditionData.decoderPosture = 'I50';
           I50_I50d_COC_block2(trial).conditionData.decoderPostureID = 1;
        end
        
    E50_E50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190729\15_externalPosture_externalDecoder_centerOutCenter_01\Earl20190729_15_externalPosture_externalDecoder_centerOutCenter_01_SI_translated.mat');
        [artifactTrials] = detectArtifact(E50_E50d_COC_block1);
        E50_E50d_COC_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        E50_E50d_COC_block1 = getDataStruct20211210(E50_E50d_COC_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E50_E50d_COC_block1,2);
        for trial = 1:numTrials
           E50_E50d_COC_block1(trial).conditionData.postureID = 2;
           E50_E50d_COC_block1(trial).conditionData.posture = 'E50';
           E50_E50d_COC_block1(trial).conditionData.decoderPosture = 'E50';
           E50_E50d_COC_block1(trial).conditionData.decoderPostureID = 2;
        end
        
    E50_I50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190729\16_externalPosture_internalDecoder_centerOutCenter\Earl20190729_16_externalPosture_internalDecoder_centerOutCenter_SI_SW_translated.mat');
    %Remove artifact
        [artifactTrials] = detectArtifact(E50_I50d_COC_block1);
        E50_I50d_COC_block1.Data(artifactTrials) = [];       
        %Get Data Struct
        E50_I50d_COC_block1 = getDataStruct20211210(E50_I50d_COC_block1.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E50_I50d_COC_block1,2);
        for trial = 1:numTrials
           E50_I50d_COC_block1(trial).conditionData.postureID = 2;
           E50_I50d_COC_block1(trial).conditionData.posture = 'E50';
           E50_I50d_COC_block1(trial).conditionData.decoderPosture = 'I50';
           E50_I50d_COC_block1(trial).conditionData.decoderPostureID = 1;
        end

    E50_E50d_COC_block2 = load('D:\Animals\Earl\2019\07\20190729\17_externalPosture_externalDecoder_centerOutCenter_02\Earl20190729_17_externalPosture_externalDecoder_centerOutCenter_02_SI_SW_translated.mat');
    %Remove artifact
        [artifactTrials] = detectArtifact(E50_E50d_COC_block2);
        E50_E50d_COC_block2.Data(artifactTrials) = [];       
        %Get Data Struct
        E50_E50d_COC_block2 = getDataStruct20211210(E50_E50d_COC_block2.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E50_E50d_COC_block2,2);
        for trial = 1:numTrials
           E50_E50d_COC_block2(trial).conditionData.postureID = 2;
           E50_E50d_COC_block2(trial).conditionData.posture = 'E50';
           E50_E50d_COC_block2(trial).conditionData.decoderPosture = 'E50';
           E50_E50d_COC_block2(trial).conditionData.decoderPostureID = 2;
        end
 
        
%% Combine Data, Keep only Successful Trials
    Data = [E50_E50d_CO_block1, I50_I50d_CO_block1, I50_I50d_COC_block1, I50_E50d_COC_block1,...
            I50_I50d_COC_block2, E50_E50d_COC_block1, E50_I50d_COC_block1,E50_E50d_COC_block2];
    
   clearvars E50_E50d_CO_block1 I50_I50d_CO_block1 I50_I50d_COC_block1 I50_E50d_COC_block1...
               I50_I50d_COC_block2 E50_E50d_COC_block1 E50_I50d_COC_block1 E50_E50d_COC_block2
    

%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    

end