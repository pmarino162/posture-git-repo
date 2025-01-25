function [Data] = loadEarlData20211004
    %Exclude Channels
    exclCh = [37:39 42:43 51 57:58 62 65:76 79:94 96 126];
    %Exclude trials
    exclTrials = [];
    
    %Load and Label Data
    F30_BC = load('D:\Animals\Earl\2021\10\20211004\06_F30_brainControl\Earl20211004_06_F30_brainControl_SI_translated.mat');
    F30_BC = F30_BC.Data;
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 4;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    
    F15_BC = load('D:\Animals\Earl\2021\10\20211004\09_F15_brainControl\Earl20211004_09_F15_brainControl_SI_translated.mat');
    F15_BC = F15_BC.Data;
    for trial = 1:size(F15_BC,2)
        F15_BC(trial).conditionData.posture = 'F15';
        F15_BC(trial).conditionData.postureID = 3;
        F15_BC(trial).conditionData.task = 'BC';
        F15_BC(trial).conditionData.taskID = 1;
    end
       
    N00_HC = load('D:\Animals\Earl\2021\10\20211004\12_N00_brainControl\Earl20211004_12_N00_brainControl_SI_translated.mat');
    N00_HC = N00_HC.Data;
    for trial = 1:size(N00_HC,2)
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'BC';
        N00_HC(trial).conditionData.taskID = 1;
    end
    
    E15_BC = load('D:\Animals\Earl\2021\10\20211004\03_E15_brainControl\Earl20211004_03_E15_brainControl_SI_translated.mat');
    E15_BC = E15_BC.Data;
    for trial = 1:size(E15_BC,2)
        E15_BC(trial).conditionData.posture = 'E15';
        E15_BC(trial).conditionData.postureID = 1;
        E15_BC(trial).conditionData.task = 'BC';
        E15_BC(trial).conditionData.taskID = 1;
    end
    
    %Get Data Struct
    Data = [F30_BC F15_BC N00_HC E15_BC];
    
    %Remove trials
    rmList = [];
    Data(rmList) = [];
    
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getKin',true,'exclCh',exclCh,'exclTrials',exclTrials); 

end