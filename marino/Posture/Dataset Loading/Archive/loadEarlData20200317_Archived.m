function [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams,postureIDs] = loadEarlData20200317(binWidth)

 %% Preprocessing Parameters 
    getMarker = true;
    getSpikes = true;
    getForce = true;
    exclCh = [3 12 16 20 31 73 85 92 94];
    
 %% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'I30';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'I15';
    postureIDs(3).ID = 3; postureIDs(3).Posture = 'N00';
    postureIDs(4).ID = 4; postureIDs(4).Posture = 'E15';
    postureIDs(5).ID = 5; postureIDs(5).Posture = 'E30';
    
%% Load, Preprocess, and Label Data    
    %3/17/2020
    date = '20200317';
%     exclCh = [];
    I30 = load('D:\Animals\Earl\2020\03\20200317\05_I30_grid\Earl20200317_05_I30_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30);
        I30.Data(artifactTrials) = [];       
        %Get Data Struct
        I30 = getDataStruct(I30.Data,'getMarker',getMarker,'getSpikes',getSpikes,'exclCh',exclCh,'binWidth',binWidth,'removeBadPhasespace',false); 
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
        I15 = getDataStruct(I15.Data,'getMarker',getMarker,'getSpikes',getSpikes,'exclCh',exclCh,'binWidth',binWidth,'removeBadPhasespace',false); 
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
        N00 = getDataStruct(N00.Data,'getMarker',getMarker,'getSpikes',getSpikes,'exclCh',exclCh,'binWidth',binWidth,'removeBadPhasespace',false); 
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
        E15 = getDataStruct(E15.Data,'getMarker',getMarker,'getSpikes',getSpikes,'exclCh',exclCh,'binWidth',binWidth,'removeBadPhasespace',false); 
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
        E30 = getDataStruct(E30.Data,'getMarker',getMarker,'getSpikes',getSpikes,'exclCh',exclCh,'binWidth',binWidth,'removeBadPhasespace',false); 
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
    noSpikeTrials = [];
    numTrials = size(Data,2);
    for trial = 1:numTrials
       allChannelSpikeBins = double(Data(trial).spikes.allChannelSpikeBins);
       sumSpikes = sum(sum(allChannelSpikeBins,1),2);
       if sumSpikes < 1
           noSpikeTrials = [noSpikeTrials,trial];
       end
    end
    Data(noSpikeTrials) = [];

    


end