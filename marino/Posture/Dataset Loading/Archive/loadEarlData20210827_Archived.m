function [Data] = loadEarlData20210827
    %Exclude Channels
    exclCh =  [35 37 42 51 57 58 62 65:68 70:76 78 81:87 89:92 126 127];
    
    %Load and Label Data
    I30_BC = load('D:\Animals\Earl\2021\08\20210827\I30_brainControl\Earl20210827_I30_brainControl_SI_translated.mat');
    I30_BC = I30_BC.Data;
    for trial = 1:size(I30_BC,2)
        I30_BC(trial).conditionData.posture = 'I30';
        I30_BC(trial).conditionData.postureID = 1;
        I30_BC(trial).conditionData.task = 'BC';
        I30_BC(trial).conditionData.taskID = 1;
    end
    
    N00_BC = load('D:\Animals\Earl\2021\08\20210827\N00_brainControl\Earl20210827_N00_brainControl_SI_translated.mat');
    N00_BC = N00_BC.Data;
    for trial = 1:size(N00_BC,2)
        N00_BC(trial).conditionData.posture = 'N00';
        N00_BC(trial).conditionData.postureID = 2;
        N00_BC(trial).conditionData.task = 'BC';
        N00_BC(trial).conditionData.taskID = 1;
    end
    
    DCO = load('D:\Animals\Earl\2021\08\20210827\DCO\Earl20210827_DCO_SI_translated.mat');
    DCO = DCO.Data;
    for trial = 1:size(DCO,2)
        DCO(trial).conditionData.posture = 'N00';
        DCO(trial).conditionData.postureID = 2;
        DCO(trial).conditionData.task = 'HC';
        DCO(trial).conditionData.taskID = 2;
    end
    
    %Get Data Struct
    Data = [I30_BC N00_BC DCO];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getForce',true,'forceSetup','shoulder_posture_bci','exclCh',exclCh); 

end