function [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314_20211210

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh =  [3 12 16 20 31 73 85 92 94];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        %BC
        checkPhasespaceSyncBC = false;
        droppedMarkerTrialsBC = 'all';
        removeBadPhasespaceBC = false;
        %HC
        checkPhasespaceSyncHC = true;
        droppedMarkerTrialsHC = 'no drop';
        removeBadPhasespaceHC = true;
        %Iso
        checkPhasespaceSyncIso = false;
        droppedMarkerTrialsIso = 'all';
        removeBadPhasespaceIso = false;
        centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = true;
        forceSetup = 'EL';
    getAlg = true;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = true;
    %BC
    checkDecodeBC = true;
    %HC
    checkDecodeHC = false;
    %Iso
    checkDecodeIso = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];
 
%% Load, Preprocess, and Label Data    
    I30_BC = load('D:\Animals\Earl\2020\03\20200314\16_internal_grid\Earl20200314_16_internal_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30_BC);
        I30_BC.Data(artifactTrials) = [];       
        %Get Data Struct
        I30_BC = getDataStruct20211210(I30_BC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(I30_BC,2);
        for trial = 1:numTrials
           I30_BC(trial).conditionData.posture = 'I30';
           I30_BC(trial).conditionData.postureID = 1;
           I30_BC(trial).conditionData.task = 'BC';
           I30_BC(trial).conditionData.taskID = 1;
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
        %Remove problem reaching trials
        [I30_HC] = cleanData20200314(I30_HC);
        %Get Data Struct
        I30_HC = getDataStruct20211210(I30_HC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncHC,'droppedMarkerTrials',droppedMarkerTrialsHC,'removeBadPhasespace',removeBadPhasespaceHC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeHC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(I30_HC,2);
        for trial = 1:numTrials
           I30_HC(trial).conditionData.posture = 'I30';
           I30_HC(trial).conditionData.postureID = 1;
           I30_HC(trial).conditionData.task = 'HC'; 
           I30_HC(trial).conditionData.taskID = 2;
        end  
    
   I30_Iso = load('D:\Animals\Earl\2020\03\20200314\18_internal_isoForce\Earl20200314_18_internal_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(I30_Iso);
        I30_Iso.Data(artifactTrials) = [];       
        %Get Data Struct
        I30_Iso = getDataStruct20211210(I30_Iso.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncIso,'droppedMarkerTrials',droppedMarkerTrialsIso,'removeBadPhasespace',removeBadPhasespaceIso,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeIso,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(I30_Iso,2);
        for trial = 1:numTrials
           I30_Iso(trial).conditionData.posture = 'I30'; 
           I30_Iso(trial).conditionData.postureID = 1; 
           I30_Iso(trial).conditionData.task = 'Iso';
           I30_Iso(trial).conditionData.taskID = 3;
        end  
        
    N00_BC = load('D:\Animals\Earl\2020\03\20200314\04_neutral_grid\Earl20200314_04_neutral_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00_BC);
        N00_BC.Data(artifactTrials) = [];
        %Get Data Struct
        N00_BC = getDataStruct20211210(N00_BC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(N00_BC,2);
        for trial = 1:numTrials
           N00_BC(trial).conditionData.posture = 'N00'; 
           N00_BC(trial).conditionData.postureID = 3;
           N00_BC(trial).conditionData.task = 'BC';
           N00_BC(trial).conditionData.taskID = 1;
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
        %Remove problem reaching trials
        [N00_HC] = cleanData20200314(N00_HC);
        %Get Data Struct
        N00_HC = getDataStruct20211210(N00_HC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncHC,'droppedMarkerTrials',droppedMarkerTrialsHC,'removeBadPhasespace',removeBadPhasespaceHC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeHC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(N00_HC,2);
        for trial = 1:numTrials
           N00_HC(trial).conditionData.posture = 'N00'; 
           N00_HC(trial).conditionData.postureID = 3; 
           N00_HC(trial).conditionData.task = 'HC';
           N00_HC(trial).conditionData.taskID = 2;
        end
   
    N00_Iso = load('D:\Animals\Earl\2020\03\20200314\06_neutral_isoForce\Earl20200314_06_neutral_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(N00_Iso);
        N00_Iso.Data(artifactTrials) = [];
        %Get Data Struct
        N00_Iso = getDataStruct20211210(N00_Iso.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncIso,'droppedMarkerTrials',droppedMarkerTrialsIso,'removeBadPhasespace',removeBadPhasespaceIso,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeIso,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(N00_Iso,2);
        for trial = 1:numTrials
           N00_Iso(trial).conditionData.posture = 'N00'; 
           N00_Iso(trial).conditionData.postureID = 3; 
           N00_Iso(trial).conditionData.task = 'Iso';
           N00_Iso(trial).conditionData.taskID = 3;
        end
        
   E30_BC = load('D:\Animals\Earl\2020\03\20200314\10_external_grid\Earl20200314_10_external_grid_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30_BC);
        E30_BC.Data(artifactTrials) = [];
        %Get Data Struct
        E30_BC = getDataStruct20211210(E30_BC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(E30_BC,2);
        for trial = 1:numTrials
           E30_BC(trial).conditionData.posture = 'E30';
           E30_BC(trial).conditionData.postureID = 5;
           E30_BC(trial).conditionData.task = 'BC'; 
           E30_BC(trial).conditionData.taskID = 1; 
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
        %Remove problem reaching trials
        [E30_HC] = cleanData20200314(E30_HC);
        %Get Data Struct
        E30_HC = getDataStruct20211210(E30_HC.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncHC,'droppedMarkerTrials',droppedMarkerTrialsHC,'removeBadPhasespace',removeBadPhasespaceHC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeHC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(E30_HC,2);
        for trial = 1:numTrials
           E30_HC(trial).conditionData.posture = 'E30';
           E30_HC(trial).conditionData.postureID = 5;
           E30_HC(trial).conditionData.task = 'HC'; 
           E30_HC(trial).conditionData.taskID = 2;
        end

    E30_Iso = load('D:\Animals\Earl\2020\03\20200314\12_external_isoForce\Earl20200314_12_external_isoForce_SI_translated.mat');
        %Remove artifact
        [artifactTrials] = detectArtifact(E30_Iso);
        E30_Iso.Data(artifactTrials) = [];
        %Get Data Struct
        E30_Iso = getDataStruct20211210(E30_Iso.Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncIso,'droppedMarkerTrials',droppedMarkerTrialsIso,'removeBadPhasespace',removeBadPhasespaceIso,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeIso,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
        numTrials = size(E30_Iso,2);
        for trial = 1:numTrials
           E30_Iso(trial).conditionData.posture = 'E30'; 
           E30_Iso(trial).conditionData.postureID = 5'; 
           E30_Iso(trial).conditionData.task = 'Iso'; 
           E30_Iso(trial).conditionData.taskID = 3;
        end
        
    Data = [I30_BC,I30_HC,I30_Iso,N00_BC,N00_HC,N00_Iso,E30_BC,E30_HC,E30_Iso];
    Data = Data([Data.trialStatus]==1);
    clearvars I30_BC I30_HC I30_Iso N00_BC N00_HC N00_Iso E30_BC E30_HC E30_Iso
    

%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);

end