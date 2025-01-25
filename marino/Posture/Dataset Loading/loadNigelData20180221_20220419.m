function [Data,postureIDs] = loadNigelData20180221_20220419

 %% Preprocessing Parameters 
    dataset = 'N20180221';
    getSpikes = true;
        getSorts = true;
        exclCh = [];
        exclZero = false; %Exclude channels that only contain sort 0
    getMarker = false;
        centerMarker = true; %Subtracts workspace center from marker pose if true
    inclStateTable = false;
    getKin = false;
    exclTrials = [];
    trialName = '';


 %% Create PostureIDs
    postureIDs = struct('ID',[],'Posture',[]);
    postureIDs(1).ID = 1; postureIDs(1).Posture = 'N00';
    postureIDs(2).ID = 2; postureIDs(2).Posture = 'I45';
    postureIDs(3).ID = 3; postureIDs(3).Posture = 'A90';

    
%% Load, Preprocess, and Label Data    
    Data = [];
    %Posture 3 (block 1)
        prevTrialNum = 0;
        for blockInd = 54306:54325
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exclTrials',exclTrials);  
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 3;
                tempData(i).conditionData.posture = 'A90';
                tempData(i).conditionData.block = 1;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
        
     %Posture 2 (block 2)
        %prevTrialNum = 0;
        for blockInd = 54331:54350
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 2;
                tempData(i).conditionData.posture = 'I45';
                tempData(i).conditionData.block = 2;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
        
      %Posture 1 (block 3)
        %prevTrialNum = 0;
        for blockInd = [54356:54369,54371:54375]
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 1;
                tempData(i).conditionData.posture = 'N00';
                tempData(i).conditionData.block = 3;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
        
      %Posture 2 (block 4)
        %prevTrialNum = 0;
        for blockInd = 54381:54400
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials);  
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 2;
                tempData(i).conditionData.posture = 'I45';
                tempData(i).conditionData.block = 4;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
        
      %Posture 3 (block 5)
        %prevTrialNum = 0;
        for blockInd = 54406:54425
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials); 
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 3;
                tempData(i).conditionData.posture = 'A90';
                tempData(i).conditionData.block = 5;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
        
      %Posture 1 (block 6)
        %prevTrialNum = 0;
        for blockInd = 54431:54525
            %Load/Preprocess Data Block
            loadStr = ['D:\Animals\Nigel\2018\02\20180221\Nigel.BC.',num2str(blockInd),'.CenterOut.mat'];
            rawData = load(loadStr);
            %[tempData] = getDataStructNigel(rawData,'getSpikes',true,'trialName','Nigel Posture BC Center Out');
            tempData = getDataStructNigel20220419(rawData,'trialName','Nigel Posture BC Center Out','getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
                'getMarker',getMarker,'centerMarker',centerMarker,'getKin',getKin,'exlcTrials',exclTrials); 
            %Update trialNum; Add Posture
            for i=1:size(tempData,2)
                tempData(i).trialNum = prevTrialNum + tempData(i).trialNum;
                tempData(i).conditionData.postureID = 1;
                tempData(i).conditionData.posture = 'N00';
                tempData(i).conditionData.block = 6;
            end
            %Add to Data Struct; Update trialNum
            Data = [Data,tempData];
            prevTrialNum = Data(end).trialNum;
        end
    


end