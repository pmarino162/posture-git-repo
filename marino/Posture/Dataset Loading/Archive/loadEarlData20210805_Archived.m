function [Data] = loadEarlData20210805
    %Exclude Channels
    exclCh = [7 12 35 40 47 58 68 70 76:78 85 87 90 126];
    
    %Load and Label Data
    N00_BC = load('D:\Animals\Earl\2021\08\20210805\17_N00_brainControl\Earl20210805_17_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 2;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    F30_BC = load('D:\Animals\Earl\2021\08\20210805\07_F30_brainControl\Earl20210805_07_F30_brainControl_SI_translated.mat');
    F30_BC = F30_BC.Data;
    for trial = 1:size(F30_BC,2)
        F30_BC(trial).conditionData.posture = 'F30';
        F30_BC(trial).conditionData.postureID = 1;
        F30_BC(trial).conditionData.task = 'BC';
        F30_BC(trial).conditionData.taskID = 1;
    end
    
    E30_BC = load('D:\Animals\Earl\2021\08\20210805\12_E30_brainControl\Earl20210805_12_E30_brainControl_SI_translated.mat');
    E30_BC = E30_BC.Data;
    for trial = 1:size(E30_BC,2)
        E30_BC(trial).conditionData.posture = 'E30';
        E30_BC(trial).conditionData.postureID = 3;
        E30_BC(trial).conditionData.task = 'BC';
        E30_BC(trial).conditionData.taskID = 1;
    end
    
    
    N00_HC = load('D:\Animals\Earl\2021\08\20210805\01_N00_activeMovements\Earl20210805_01_N00_activeMovements_SI_translated.mat');
    N00_HC = N00_HC.Data;
    %Keep only trials with correct phasespace-vr transformation
    keepList = [];
    for trial = 1:size(N00_HC,2)
        if N00_HC(trial).Definitions.vrPhasespaceTransformation.translation == [0 390 -680]
            keepList = [keepList,trial];
        end
        N00_HC(trial).conditionData.posture = 'N00';
        N00_HC(trial).conditionData.postureID = 2;
        N00_HC(trial).conditionData.task = 'HC';
        N00_HC(trial).conditionData.taskID = 2;
    end
    N00_HC = N00_HC(keepList);
    
    E30_HC = load('D:\Animals\Earl\2021\08\20210805\02_E30_activeMovements\Earl20210805_02_E30_activeMovements_SI_translated.mat');
    E30_HC = E30_HC.Data;
    %Keep only trials with correct phasespace-vr transformation
    keepList = [];
    for trial = 1:size(E30_HC,2)
        if E30_HC(trial).Definitions.vrPhasespaceTransformation.translation == [-20 520 -680]
            keepList = [keepList,trial];
        end
        E30_HC(trial).conditionData.posture = 'E30';
        E30_HC(trial).conditionData.postureID = 3;
        E30_HC(trial).conditionData.task = 'HC';
        E30_HC(trial).conditionData.taskID = 2;
    end
    E30_HC = E30_HC(keepList);

    %Get Data Struct
    Data = [N00_BC F30_BC E30_BC, N00_HC, E30_HC];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'exclCh',exclCh); 

end