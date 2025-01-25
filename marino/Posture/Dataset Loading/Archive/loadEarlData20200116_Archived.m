function [Data] = loadEarlData20200116(binWidth)

    %Exclude Channels
    exclCh = [3 12 16 20 31 62 73 94];
    %Exclude trials
    exclTrials = [];
    
    %Load and Label Data
    I30_Iso = load('D:\Animals\Earl\2020\01\20200116\06_isoForceI30\Earl20200116_06_isoForceI30_SI_SW_translated.mat');
    I30_Iso = I30_Iso.Data;
    for trial = 1:size(I30_Iso,2)
        I30_Iso(trial).conditionData.posture = 'I30';
        I30_Iso(trial).conditionData.postureID = 1;
        I30_Iso(trial).conditionData.task = 'Iso';
        I30_Iso(trial).conditionData.taskID = 3;
    end
    
    I15_Iso =  load('D:\Animals\Earl\2020\01\20200116\05_isoForceI15\Earl20200116_05_isoForceI15_SI_SW_translated.mat');
    I15_Iso = I15_Iso.Data;
    for trial = 1:size(I15_Iso,2)
        I15_Iso(trial).conditionData.posture = 'I15';
        I15_Iso(trial).conditionData.postureID = 2;
        I15_Iso(trial).conditionData.task = 'Iso';
        I15_Iso(trial).conditionData.taskID = 3;
    end
    
    
    N00_Iso = load('D:\Animals\Earl\2020\01\20200116\07_isoForceN\Earl20200116_07_isoForceN_SI_SW_translated.mat');
    N00_Iso = N00_Iso.Data;
    for trial = 1:size(N00_Iso,2)
        N00_Iso(trial).conditionData.posture = 'N00';
        N00_Iso(trial).conditionData.postureID = 3;
        N00_Iso(trial).conditionData.task = 'Iso';
        N00_Iso(trial).conditionData.taskID = 3;
    end
    
    E15_Iso =  load('D:\Animals\Earl\2020\01\20200116\03_isoForceE15\Earl20200116_03_isoForceE15_SI_SW_translated.mat');
    E15_Iso = E15_Iso.Data;
    for trial = 1:size(E15_Iso,2)
        E15_Iso(trial).conditionData.posture = 'E15';
        E15_Iso(trial).conditionData.postureID = 4;
        E15_Iso(trial).conditionData.task = 'Iso';
        E15_Iso(trial).conditionData.taskID = 3;
    end
        
    E30_Iso =  load('D:\Animals\Earl\2020\01\20200116\04_isoForceE30\Earl20200116_04_isoForceE30_SI_SW_translated.mat');
    E30_Iso = E30_Iso.Data;
    for trial = 1:size(E30_Iso,2)
        E30_Iso(trial).conditionData.posture = 'E30';
        E30_Iso(trial).conditionData.postureID = 5;
        E30_Iso(trial).conditionData.task = 'Iso';
        E30_Iso(trial).conditionData.taskID = 3;
    end
    
    %Get Data Struct, Keep only successful trials
    Data = [I30_Iso I15_Iso N00_Iso E15_Iso E30_Iso];
    Data = getDataStruct(Data,'getMarker',true,'getSpikes',true,'getKin',true,'exclCh',exclCh,'binWidth',binWidth,'exclTrials',exclTrials); 
    Data = Data([Data.trialStatus]==1);
end