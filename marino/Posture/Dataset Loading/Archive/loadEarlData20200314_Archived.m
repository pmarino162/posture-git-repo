function [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadData20200314

%% Load, Preprocess, and Label Data    
    %3/14/2020
    date = '20200314';
    exclCh =  [3 12 16 20 31 73 85 92 94];

    I30_BC = load('D:\Animals\Earl\2020\03\20200314\16_internal_grid\Earl20200314_16_internal_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30_BC);
        I30_BC.Data(artifactTrials) = [];       
        %Get Data Struct
        I30_BC = getDataStruct(I30_BC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(I30_BC,2);
        for trial = 1:numTrials
           I30_BC(trial).postureData.posture = 'I30';
           I30_BC(trial).postureData.postureID = 1;
           I30_BC(trial).taskData.task = 'BC';
           I30_BC(trial).taskData.taskID = 1;
           I30_BC(trial).targetData.targetID = I30_BC(trial).targetData.target1ID;
        end
        %Load GPFA and Decoder Params
        I30GPFAParams = load('D:\Animals\Earl\2020\03\20200314\analysis\mat_results\run015\gpfa_xDim10.mat');
        I30DecoderParams = I30_BC(1).Decoder.Parameters;        
   
    I30_HC = load('D:\Animals\Earl\2020\03\20200314\17_internal_forceBar_CO\Earl20200314_17_internal_forceBar_CO_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30_HC);
        I30_HC.Data(artifactTrials) = [];      
        %Rename trials
        numTrials = size(I30_HC.Data,2);
        for trial = 1:numTrials
            I30_HC.Data(trial).Overview.trialName = 'HC_CenterOut_ForceBar_20200314';
        end
        %Get Data Struct
        I30_HC = getDataStruct(I30_HC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(I30_HC,2);
        for trial = 1:numTrials
           I30_HC(trial).postureData.posture = 'I30';
           I30_HC(trial).postureData.postureID = 1;
           I30_HC(trial).taskData.task = 'HC'; 
           I30_HC(trial).taskData.taskID = 2;
        end  
    
   I30_Iso = load('D:\Animals\Earl\2020\03\20200314\18_internal_isoForce\Earl20200314_18_internal_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30_Iso);
        I30_Iso.Data(artifactTrials) = [];       
        %Get Data Struct
        I30_Iso = getDataStruct(I30_Iso.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(I30_Iso,2);
        for trial = 1:numTrials
           I30_Iso(trial).postureData.posture = 'I30'; 
           I30_Iso(trial).postureData.postureID = 1; 
           I30_Iso(trial).taskData.task = 'Iso';
           I30_Iso(trial).taskData.taskID = 3;
        end  
        
    N00_BC = load('D:\Animals\Earl\2020\03\20200314\04_neutral_grid\Earl20200314_04_neutral_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00_BC);
        N00_BC.Data(artifactTrials) = [];
        %Get Data Struct
        N00_BC = getDataStruct(N00_BC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(N00_BC,2);
        for trial = 1:numTrials
           N00_BC(trial).postureData.posture = 'N00'; 
           N00_BC(trial).postureData.postureID = 3;
           N00_BC(trial).taskData.task = 'BC';
           N00_BC(trial).taskData.taskID = 1;
           N00_BC(trial).targetData.targetID = N00_BC(trial).targetData.target1ID;
        end
        %Load GPFA and Decoder Params
        N00GPFAParams = load('D:\Animals\Earl\2020\03\20200314\analysis\mat_results\run005\gpfa_xDim10.mat');
        N00DecoderParams = N00_BC(1).Decoder.Parameters;          
        
    N00_HC = load('D:\Animals\Earl\2020\03\20200314\05_neutral_forceBar_CO\Earl20200314_05_neutral_forceBar_CO_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00_HC);
        N00_HC.Data(artifactTrials) = [];
        %Rename trials
        numTrials = size(N00_HC.Data,2);
        for trial = 1:numTrials
            N00_HC.Data(trial).Overview.trialName = 'HC_CenterOut_ForceBar_20200314';
        end
        %Get Data Struct
        N00_HC = getDataStruct(N00_HC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(N00_HC,2);
        for trial = 1:numTrials
           N00_HC(trial).postureData.posture = 'N00'; 
           N00_HC(trial).postureData.postureID = 3; 
           N00_HC(trial).taskData.task = 'HC';
           N00_HC(trial).taskData.taskID = 2;
        end
   
    N00_Iso = load('D:\Animals\Earl\2020\03\20200314\06_neutral_isoForce\Earl20200314_06_neutral_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00_Iso);
        N00_Iso.Data(artifactTrials) = [];
        %Get Data Struct
        N00_Iso = getDataStruct(N00_Iso.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(N00_Iso,2);
        for trial = 1:numTrials
           N00_Iso(trial).postureData.posture = 'N00'; 
           N00_Iso(trial).postureData.postureID = 5; 
           N00_Iso(trial).taskData.task = 'Iso';
           N00_Iso(trial).taskData.taskID = 3;
        end
        
   E30_BC = load('D:\Animals\Earl\2020\03\20200314\10_external_grid\Earl20200314_10_external_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30_BC);
        E30_BC.Data(artifactTrials) = [];
        %Get Data Struct
        E30_BC = getDataStruct(E30_BC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(E30_BC,2);
        for trial = 1:numTrials
           E30_BC(trial).postureData.posture = 'E30';
           E30_BC(trial).postureData.postureID = 5;
           E30_BC(trial).taskData.task = 'BC'; 
           E30_BC(trial).taskData.taskID = 1; 
           E30_BC(trial).targetData.targetID = E30_BC(trial).targetData.target1ID;
        end
        %Load GPFA and Decoder Params
        E30GPFAParams = load('D:\Animals\Earl\2020\03\20200314\analysis\mat_results\run010\gpfa_xDim10.mat');
        E30DecoderParams = E30_BC(1).Decoder.Parameters;
        
    E30_HC = load('D:\Animals\Earl\2020\03\20200314\11_external_forceBar_CO\Earl20200314_11_external_forceBar_CO_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30_HC);
        E30_HC.Data(artifactTrials) = [];
        %Rename trials
        numTrials = size(E30_HC.Data,2);
        for trial = 1:numTrials
            E30_HC.Data(trial).Overview.trialName = 'HC_CenterOut_ForceBar_20200314';
        end
        %Get Data Struct
        E30_HC = getDataStruct(E30_HC.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(E30_HC,2);
        for trial = 1:numTrials
           E30_HC(trial).postureData.posture = 'E30';
           E30_HC(trial).postureData.postureID = 5;
           E30_HC(trial).taskData.task = 'HC'; 
           E30_HC(trial).taskData.taskID = 2;
        end

    E30_Iso = load('D:\Animals\Earl\2020\03\20200314\12_external_isoForce\Earl20200314_12_external_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30_Iso);
        E30_Iso.Data(artifactTrials) = [];
        %Get Data Struct
        E30_Iso = getDataStruct(E30_Iso.Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 
        numTrials = size(E30_Iso,2);
        for trial = 1:numTrials
           E30_Iso(trial).postureData.posture = 'E30'; 
           E30_Iso(trial).postureData.postureID = 5'; 
           E30_Iso(trial).taskData.task = 'Iso'; 
           E30_Iso(trial).taskData.taskID = 3;
        end
        
    Data = [I30_BC,I30_HC,I30_Iso,N00_BC,N00_HC,N00_Iso,E30_BC,E30_HC,E30_Iso];
    Data = Data([Data.trialStatus]==1);
    clearvars I30_BC I30_HC I30_Iso N00_BC N00_HC N00_Iso E30_BC E30_HC E30_Iso
    

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


% %% Label Conditions
%     numTrials = size(Data,2);
%     for trial = 1:numTrials
%        Data(trial).Condition = [Data(trial).Posture,'_',Data(trial).Task,'_T',num2str(Data(trial).targetData.targetID)]; 
%     end

    
% %% Get minimum number of successful trials for any condition
%     conditionList = unique({Data.Condition});
%     minNumber = size(Data,2);
%     maxNumber = 0;
%     for condition = conditionList
%        tempData = Data(strcmpi({Data.Condition},condition))
%        trialStatus = [tempData.trialStatus];
%        numSuc = sum(trialStatus);
%        if numSuc < minNumber
%            minNumber = numSuc;
%        end
%        if numSuc > maxNumber
%            maxNumber = numSuc;
%        end
%     end
%     
%     
% %% Keep only minimum number
%     minSucData = [];
%     for condition = conditionList
%         tempData = Data(strcmpi({Data.Condition},condition));
%         minSucData = [minSucData,tempData(1:minNumber)];
%     end
%     Data = minSucData;
%     clearvars minSucData


end