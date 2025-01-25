function [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams,postureIDs] = loadEarlData20200317_20211210

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh =  [3 12 16 20 31 73 85 92 94];
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
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'I30';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'I15';
    postureIDs(3).ID = 3; postureIDs(3).Posture = 'N00';
    postureIDs(4).ID = 4; postureIDs(4).Posture = 'E15';
    postureIDs(5).ID = 5; postureIDs(5).Posture = 'E30';
    
%% Load, Preprocess, and Label Data    
    I30 = load('D:\Animals\Earl\2020\03\20200317\05_I30_grid\Earl20200317_05_I30_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30);
        I30.Data(artifactTrials) = [];       
        %Get Data Struct
        I30 = getDataStruct20211210(I30.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I30,2);
        for trial = 1:numTrials
           I30(trial).conditionData.postureID = 1;
           I30(trial).conditionData.posture = 'I30';
        end
        %Load GPFA and Decoder Params
        I30GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run005\gpfa_xDim10.mat');
        I30DecoderParams = I30(1).Decoder.Parameters;
        
    I15 = load('D:\Animals\Earl\2020\03\20200317\09_I15_grid\Earl20200317_09_I15_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I15);
        I15.Data(artifactTrials) = [];
        %Get Data Struct
        I15 = getDataStruct20211210(I15.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(I15,2);
        for trial = 1:numTrials
           I15(trial).conditionData.postureID = 2;
           I15(trial).conditionData.posture = 'I15';
        end
        %Load GPFA and Decoder Params
        I15GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run010\gpfa_xDim10.mat');
        I15DecoderParams = I15(1).Decoder.Parameters;
        
    N00 = load('D:\Animals\Earl\2020\03\20200317\13_N_grid\Earl20200317_13_N_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00);
        N00.Data(artifactTrials) = [];
        %Get Data Struct
        N00 = getDataStruct20211210(N00.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(N00,2);
        for trial = 1:numTrials
           N00(trial).conditionData.postureID = 3;
           N00(trial).conditionData.posture = 'N00';
        end
        %Load GPFA and Decoder Params
        N00GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run015\gpfa_xDim10.mat');
        N00DecoderParams = N00(1).Decoder.Parameters;
        
    E15 = load('D:\Animals\Earl\2020\03\20200317\17_E15_grid\Earl20200317_17_E15_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E15);
        E15.Data(artifactTrials) = [];
        %Get Data Struct
        E15 = getDataStruct20211210(E15.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E15,2);
        for trial = 1:numTrials
           E15(trial).conditionData.postureID = 4;
           E15(trial).conditionData.posture = 'E15';
        end
        %Load GPFA and Decoder Params
        E15GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run020\gpfa_xDim10.mat');
        E15DecoderParams = E15(1).Decoder.Parameters;
        
   E30 = load('D:\Animals\Earl\2020\03\20200317\21_E30_grid\Earl20200317_21_E30_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30);
        E30.Data(artifactTrials) = [];
        %Get Data Struct
        E30 = getDataStruct20211210(E30.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);      
        numTrials = size(E30,2);
        for trial = 1:numTrials
           E30(trial).conditionData.postureID = 5;
           E30(trial).conditionData.posture = 'E30';
        end
        %Load GPFA and Decoder Params
        E30GPFAParams = load('D:\Animals\Earl\2020\03\20200317\analysis\mat_results\run025\gpfa_xDim10.mat');
        E30DecoderParams = E30(1).Decoder.Parameters;
    
%% Combine Data, Keep only Successful Trials
    Data = [I30,I15,N00,E15,E30];
    Data = Data([Data.trialStatus]==1);
    clearvars I30 I15 N00 E15 E30
    

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