function [Data] = loadEarlData20210901
    %Exclude Channels
    exclCh = [35 39 43 51 56:58 65:78 81:87 89:93 96];
    %Exclude trials
    exclTrials = [];
    
    %Load and Label Data
    F30_BC = load('D:\Animals\Earl\2021\09\20210901\05_F30_brainControl\Earl20210901_05_F30_brainControl_SI_translated.mat');
    F30_BC = F30_BC.Data;
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 4;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    
    eE30_BC = load('D:\Animals\Earl\2021\09\20210901\08_eE30_brainControl\Earl20210901_08_eE30_brainControl_SI_translated.mat');
    eE30_BC = eE30_BC.Data;
    for trial = 1:size(eE30_BC,2)
        eE30_BC(trial).conditionData.posture = 'eE30';
        eE30_BC(trial).conditionData.postureID = 5;
        eE30_BC(trial).conditionData.task = 'BC';
        eE30_BC(trial).conditionData.taskID = 1;
    end
    
    I30_BC = load('D:\Animals\Earl\2021\09\20210901\14_I30_brainControl\Earl20210901_14_I30_brainControl_SI_translated.mat');
    I30_BC = I30_BC.Data;
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 1;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    E30_BC = load('D:\Animals\Earl\2021\09\20210901\11_E30_brainControl\Earl20210901_11_E30_brainControl_SI_translated.mat');
    E30_BC = E30_BC.Data;
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 3;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
       
    N00_HC = load('D:\Animals\Earl\2021\09\20210901\02_DCO\Earl20210901_02_DCO_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
    %Get Data Struct
    Data = [F30_BC eE30_BC I30_BC E30_BC N00_HC];
    
    %Remove alignment trials
    rmList = [];
    for i = 1:size(Data,2)
        if strcmpi(Data(i).Parameters.TrialTargets.names{1,3},'chair alignment target')
            rmList = [rmList,i];
        end
    end
    Data(rmList) = [];
    
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getKin',true,'exclCh',exclCh,'exclTrials',exclTrials); 

end