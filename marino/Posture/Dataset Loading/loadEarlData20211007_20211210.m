function [Data] = loadEarlData20211007_20211210()

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh = [37:39, 42, 43, 51, 53, 57:58, 62, 65:76, 79:94, 96, 126];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = false;
        centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = false;
        forceSetup = '';
    getAlg = false;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = false;
    checkDecode = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];

%% Load, Preprocess, and Label Data 
    N00_Iso = load('D:\Animals\Earl\2021\10\20211007\N00_isoForce\Earl20211007_N00_isoForce_SI_translated.mat');
        %Remove artifact
        %[artifactTrials] = detectArtifact(N00_Iso);
        %N00_Iso.Data(artifactTrials) = [];
        N00_Iso = N00_Iso.Data;
        for trial = 1:size(N00_Iso,2)
            N00_Iso(trial).conditionData.posture = 'N00';
            N00_Iso(trial).conditionData.postureID = 1;
            N00_Iso(trial).conditionData.task = 'Iso';
            N00_Iso(trial).conditionData.taskID = 3;
        end

    I15_Iso = load('D:\Animals\Earl\2021\10\20211007\I15_isoForce\Earl20211007_I15_isoForce_SI_translated.mat');
        %Remove artifact
        %[artifactTrials] = detectArtifact(I15_Iso);
        %I15_Iso.Data(artifactTrials) = [];
        I15_Iso = I15_Iso.Data;
        for trial = 1:size(I15_Iso,2)
            I15_Iso(trial).conditionData.posture = 'I15';
            I15_Iso(trial).conditionData.postureID = 2;
            I15_Iso(trial).conditionData.task = 'Iso';
            I15_Iso(trial).conditionData.taskID = 3;
        end    
                        
    I30_Iso = load('D:\Animals\Earl\2021\10\20211007\I30_isoForce\Earl20211007_I30_isoForce_SI_translated.mat');
        %Remove artifact
        %[artifactTrials] = detectArtifact(I30_Iso);
        %I30_Iso.Data(artifactTrials) = [];
        I30_Iso = I30_Iso.Data;
        for trial = 1:size(I30_Iso,2)
            I30_Iso(trial).conditionData.posture = 'I30';
            I30_Iso(trial).conditionData.postureID = 3;
            I30_Iso(trial).conditionData.task = 'Iso';
            I30_Iso(trial).conditionData.taskID = 3;
        end
    
        I45_Iso = load('D:\Animals\Earl\2021\10\20211007\I45_isoForce\Earl20211007_I45_isoForce_SI_translated.mat');
        %Remove artifact
        %[artifactTrials] = detectArtifact(I45_Iso);
        %I45_Iso.Data(artifactTrials) = [];
        I45_Iso = I45_Iso.Data;
        for trial = 1:size(I45_Iso,2)
            I45_Iso(trial).conditionData.posture = 'I45';
            I45_Iso(trial).conditionData.postureID = 4;
            I45_Iso(trial).conditionData.task = 'Iso';
            I45_Iso(trial).conditionData.taskID = 3;
        end

    
%% Get Data Struct, Keep only successful trials
    Data = [N00_Iso I15_Iso I30_Iso I45_Iso];

    Data = getDataStruct20211210(Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);  
    
    Data = Data([Data.trialStatus]==1);
    
%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
end