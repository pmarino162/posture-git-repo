function [Data] = loadEarlData20210828
    %Exclude Channels
    exclCh = [35 37 42 51 57 58 62 65:68 70:76 78 81:87 89:92 126 127]; 
    %Exclude trials
    exclTrials = [];
    
    %Load and Label Data
    N00_BC = load('D:\Animals\Earl\2021\08\20210828\N00_brainControl\Earl20210828_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 2;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    I30_BC = load('D:\Animals\Earl\2021\08\20210828\I30_brainControl\Earl20210828_I30_brainControl_SI_translated.mat');
    I30_BC = I30_BC.Data;
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 1;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    E30_BC = load('D:\Animals\Earl\2021\08\20210828\E30_brainControl\Earl20210828_E30_brainControl_SI_translated.mat');
    E30_BC = E30_BC.Data;
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 3;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
       
    N00_HC = load('D:\Animals\Earl\2021\08\20210828\DCO\Earl20210828_DCO_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
    %Get Data Struct
    Data = [I30_BC N00_BC E30_BC N00_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getKin',true,'exclCh',exclCh,'exclTrials',exclTrials); 

end