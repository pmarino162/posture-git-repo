function [Data] = loadEarlData20210825
    %Exclude Channels
    exclCh = [7 35 42 43 47 51 57 58 62 65:78 81:87 89:92 118 126 127];
    
    %Load and Label Data
    BC = load('D:\Animals\Earl\2021\08\20210825\N00_brainControl\Earl20210825_N00_brainControl_SI_translated.mat');
    BC = BC.Data;
    for trial = 1:size(BC,2)
        BC(trial).conditionData.task = 'BC';
        BC(trial).conditionData.taskID = 1;
    end
    
    DCO = load('D:\Animals\Earl\2021\08\20210825\DCO\Earl20210825_DCO_SI_translated.mat');
    DCO = DCO.Data;
    for trial = 1:size(DCO,2)
        DCO(trial).conditionData.task = 'HC';
        DCO(trial).conditionData.taskID = 2;
    end
    
    %Get Data Struct
    Data = [BC DCO];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getForce',true,'forceSetup','shoulder_posture_bci','exclCh',exclCh); 

end