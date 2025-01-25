function [Data,postureIDs] = loadEarlData20190728

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh =  [];%[3 12 16 20 31 73 85 92 94];
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
    I50_I50d_CO_block1 = load('D:\Animals\Earl\2019\07\20190728\05_internal_centerOut_01\Earl20190728_05_internal_centerOut_01_SI_SW_translated.mat');
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
    
    E50_E50d_CO_block1 = load('D:\Animals\Earl\2019\07\20190728\10_external_centerOut\Earl20190728_10_external_centerOut_SI_SW_translated.mat');
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


    E50_I50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190728\12_externalPosture_internalDecoder_centerOutCenter\Earl20190728_12_externalPosture_internalDecoder_centerOutCenter_SI_SW_translated.mat');
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

    E50_E50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190728\13_externalPosture_externalDecoder_centerOutCenter\Earl20190728_13_externalPosture_externalDecoder_centerOutCenter_SI_SW_translated.mat');
    %Remove artifact
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


    I50_E50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190728\14_internalPosture_externalDecoder_centerOutCenter_01\Earl20190728_14_internalPosture_externalDecoder_centerOutCenter_SI_SW_translated.mat');    
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

        I50_I50d_COC_block1 = load('D:\Animals\Earl\2019\07\20190728\15_internalPosture_internalDecoder_centerOutCenter_01\Earl20190728_15_internalPosture_internalDecoder_centerOutCenter_SI_SW_translated.mat');
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
           I50_I50d_COC_block1(trial).conditionData.decoderPosture = 'E50';
           I50_I50d_COC_block1(trial).conditionData.decoderPostureID = 2;
        end
        
        I50_E50d_COC_block2 = load('D:\Animals\Earl\2019\07\20190728\16_internalPosture_externalDecoder_centerOutCenter_02\Earl20190728_16_internalPosture_externalDecoder_centerOutCenter_02_SI_SW_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I50_E50d_COC_block2);
        I50_E50d_COC_block2.Data(artifactTrials) = [];       
        %Get Data Struct
        I50_E50d_COC_block2 = getDataStruct20211210(I50_E50d_COC_block2.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I50_E50d_COC_block2,2);
        for trial = 1:numTrials
           I50_E50d_COC_block2(trial).conditionData.postureID = 1;
           I50_E50d_COC_block2(trial).conditionData.posture = 'I50';
           I50_E50d_COC_block2(trial).conditionData.decoderPosture = 'E50';
           I50_E50d_COC_block2(trial).conditionData.decoderPostureID = 2;
        end
        
        I50_I50d_COC_block2 = load('D:\Animals\Earl\2019\07\20190728\17_internalPosture_internalDecoder_centerOutCenter_02\Earl20190728_17_internalPosture_internalDecoder_centerOutCenter_02_SI_SW_translated.mat');
        %Remove artifact
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
           I50_I50d_COC_block2(trial).conditionData.decoderPosture = 'E50';
           I50_I50d_COC_block2(trial).conditionData.decoderPostureID = 2;
        end
 
%% Combine Data, Keep only Successful Trials
    Data = [I50_I50d_CO_block1 E50_E50d_CO_block1 E50_I50d_COC_block1 E50_E50d_COC_block1...
        I50_E50d_COC_block1 I50_I50d_COC_block1 I50_E50d_COC_block2 I50_I50d_COC_block2];
    
    Data = Data([Data.trialStatus]==1);
    clearvars I50_I50d_CO_block1 E50_E50d_CO_block1 E50_I50d_COC_block1 E50_E50d_COC_block1...
        I50_E50d_COC_block1 I50_I50d_COC_block1 I50_E50d_COC_block2 I50_I50d_COC_block2
    

%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    
% %% Exclude trials without step1AcqTime or with step1AcqTime > 750ms
%     rmTrials = [];
%     for i = 1:size(Data,2)
%         if isempty(Data(i).kinData)
%             rmTrials = [rmTrials,i];
%         elseif isempty(Data(i).kinData.step1AcqTime)
%             rmTrials = [rmTrials,i];
%         elseif Data(i).kinData.step1AcqTime > 750
%             rmTrials = [rmTrials,i];
%         end
%     end
%     
%     Data(rmTrials) = [];

end