function [Data] = loadEarlData20210720
    %Exclude Channels
    exclCh =  [7 55 67 68 70 76 77 82 87 116];
    
    %Load and Label Data
    E30_BC = load('D:\Animals\Earl\2021\07\20210720\10_E30_brainControl\Earl20210720_10_E30_brainControl_SI_translated.mat');
    E30_BC = E30_BC.Data;
    %Remove trials when Earl's arm was out of the posture device 
    E30_BC(1:22) = [];
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 1;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
    
    I30_BC = load('D:\Animals\Earl\2021\07\20210720\07_I30_brainControl\Earl20210720_07_I30_brainControl_SI_translated.mat');
    I30_BC = I30_BC.Data;
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 3;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    I60_BC = load('D:\Animals\Earl\2021\07\20210720\04_I60_brainControl\Earl20210720_04_I60_brainControl_SI_translated.mat');
    I60_BC = I60_BC.Data;
    for trial = 1:size(I60_BC,2)
        I60_BC(trial).conditionData.posture = 'I60';
        I60_BC(trial).conditionData.postureID = 4;
        I60_BC(trial).conditionData.task = 'BC';
        I60_BC(trial).conditionData.taskID = 1;
    end

    N00_HC = load('D:\Animals\Earl\2021\07\20210720\01_N00_activeMovements\Earl20210720_01_N00_activeMovements_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    
    
    %Get Data Struct
    Data = [E30_BC I30_BC I60_BC N00_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 

end