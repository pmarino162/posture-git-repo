function [Data] = loadEarlData20210901_20211210

 %% Preprocessing Parameters 
    getSpikes = true;
        getSorts = false;
        exclCh = [35 39 43 51 56:58 65:78 81:87 89:93 96];
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
    F30_BC = load('D:\Animals\Earl\2021\09\20210901\05_F30_brainControl\Earl20210901_05_F30_brainControl_SI_translated.mat');
    F30_BC = removeAlignmentTrials(F30_BC.Data);
    F30_BC = getDataStruct20211210(F30_BC,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exclTrials',exclTrials);        
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 4;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    
    eE30_BC = load('D:\Animals\Earl\2021\09\20210901\08_eE30_brainControl\Earl20210901_08_eE30_brainControl_SI_translated.mat');
    eE30_BC = removeAlignmentTrials(eE30_BC.Data);
    eE30_BC = getDataStruct20211210(eE30_BC,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);        
    for trial = 1:size(eE30_BC,2)
        eE30_BC(trial).conditionData.posture = 'eE30';
        eE30_BC(trial).conditionData.postureID = 5;
        eE30_BC(trial).conditionData.task = 'BC';
        eE30_BC(trial).conditionData.taskID = 1;
    end
    
    I30_BC = load('D:\Animals\Earl\2021\09\20210901\14_I30_brainControl\Earl20210901_14_I30_brainControl_SI_translated.mat');
    I30_BC = removeAlignmentTrials(I30_BC.Data);
    I30_BC = getDataStruct20211210(I30_BC,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);        
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 1;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    E30_BC = load('D:\Animals\Earl\2021\09\20210901\11_E30_brainControl\Earl20210901_11_E30_brainControl_SI_translated.mat');
    E30_BC = removeAlignmentTrials(E30_BC.Data);
    E30_BC = getDataStruct20211210(E30_BC,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncBC,'droppedMarkerTrials',droppedMarkerTrialsBC,'removeBadPhasespace',removeBadPhasespaceBC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeBC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);        
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 3;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
       
    N00_HC = load('D:\Animals\Earl\2021\09\20210901\02_DCO\Earl20210901_02_DCO_SI_translated.mat');
    N00_HC = removeAlignmentTrials(N00_HC.Data);
    N00_HC = getDataStruct20211210(N00_HC,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
            'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSyncHC,'droppedMarkerTrials',droppedMarkerTrialsHC,'removeBadPhasespace',removeBadPhasespaceHC,'centerMarker',centerMarker,...
            'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
            'checkDecode',checkDecodeHC,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);  
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
%% Combine data; keep only successful
    Data = [F30_BC eE30_BC I30_BC E30_BC N00_HC];
    Data = Data([Data.trialStatus]==1);
    
%% Local function for removing alignment trials
    function rawData = removeAlignmentTrials(rawData)
        rmList = [];
        for i = 1:size(rawData,2)
            if strcmpi(rawData(i).Parameters.TrialTargets.names{1,3},'chair alignment target')
                rmList = [rmList,i];
            end
        end
        rawData(rmList) = [];
    end
        

%% Remove "no-spike" trials
    Data = removeNoSpikeTrials(Data);
    
end