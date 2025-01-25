function [Data] = loadEarlData20210809
    %Exclude Channels
    exclCh = [7 12 43 44 47 50 56 58 60 67 68 70 74 76:78 81 83:85 87 90 126]; 
    %Exclude trials
    exclTrials = [1013,1030];
    
    %Load and Label Data
    N00_BC = load('D:\Animals\Earl\2021\08\20210809\14_N00_brainControl\Earl20210809_14_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 2;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    F30_BC = load('D:\Animals\Earl\2021\08\20210809\04_F30_brainControl\Earl20210809_04_F30_brainControl_SI_translated.mat');
    F30_BC = F30_BC.Data;
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 1;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    
    E30_BC = load('D:\Animals\Earl\2021\08\20210809\09_E30_brainControl\Earl20210809_09_E30_brainControl_SI_translated.mat');
    E30_BC = E30_BC.Data;
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 3;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
       
    N00_HC = load('D:\Animals\Earl\2021\08\20210809\01_N00_activeMovements\Earl20210809_01_N00_activeMovements_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
    %Get Data Struct
    Data = [F30_BC N00_BC E30_BC N00_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh,'exclTrials',exclTrials); 

end