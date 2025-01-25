function [Data] = loadEarlData20210804
    %Exclude Channels
    exclCh = [7 12 35 40 67 68 70 74 76:79 82 87 90 101 126 128];
    
    %Load and Label Data
    N00_BC = load('D:\Animals\Earl\2021\08\20210804\07_N00_brainControl\Earl20210804_07_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 1;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    F30_BC = load('D:\Animals\Earl\2021\08\20210804\04_F30_brainControl\Earl20210804_04_F30_brainControl_SI_translated.mat');
    F30_BC = F30_BC.Data;
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 2;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    

    N00_HC = load('D:\Animals\Earl\2021\08\20210804\01_N00_activeMovements\Earl20210804_01_N00_activeMovements_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 1;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
    %Get Data Struct
    Data = [N00_BC F30_BC N00_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 

end