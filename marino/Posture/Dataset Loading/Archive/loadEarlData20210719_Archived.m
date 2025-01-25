function [Data] = loadEarlData20210719
    %Exclude Channels
    exclCh = [7 67 68 82 87 116];
    
    %Load and Label Data
    N00_BC = load('D:\Animals\Earl\2021\07\20210719\04_N00_brainControl\Earl20210719_04_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 1;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    
    I30_BC = load('D:\Animals\Earl\2021\07\20210719\07_I30_brainControl\Earl20210719_07_I30_brainControl_SI_translated.mat');
    I30_BC = I30_BC.Data;
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 2;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    I60_BC = load('D:\Animals\Earl\2021\07\20210719\10_I60_brainControl\Earl20210719_10_I60_brainControl_SI_translated.mat');
    I60_BC = I60_BC.Data;
    for trial = 1:size(I60_BC,2)
        I60_BC(trial).conditionData.posture = 'I60';
        I60_BC(trial).conditionData.postureID = 3;
        I60_BC(trial).conditionData.task = 'BC';
        I60_BC(trial).conditionData.taskID = 1;
    end

    N00_HC = load('D:\Animals\Earl\2021\07\20210719\01_N00_activeMovements\Earl20210719_01_N00_activeMovements_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 1;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    I30_HC = load('D:\Animals\Earl\2021\07\20210719\11_I30_activeMovements\Earl20210719_11_I30_activeMovements_translated.mat');
    I30_HC = I30_HC.Data;
    for trial = 1:size(I30_HC,2)
        I30_HC(trial).conditionData.posture = 'I30';
        I30_HC(trial).conditionData.postureID = 2;
        I30_HC(trial).conditionData.task = 'HC';
        I30_HC(trial).conditionData.taskID = 2;
    end
    
    %Get Data Struct
    Data = [N00_BC I30_BC I60_BC N00_HC I30_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 

end