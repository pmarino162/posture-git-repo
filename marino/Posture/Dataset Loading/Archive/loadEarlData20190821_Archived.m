function [Data,GPFAParams,DecoderParams] = loadEarlData20190821

%% Load, Preprocess, and Label Data    
    %8/21/2019
    date = '20190821';

    BC = load('D:\Animals\Earl\2020\08\20190821\04_centerOutCenter_BC\Earl20190821_04_centerOutCenter_BC_SI_translated.mat');
        %Get Data Struct
        BC = getDataStruct(BC.Data,'getMarker',true,'getSpikes',true); 
        numTrials = size(BC,2);
        for trial = 1:numTrials
           BC(trial).Task = 'BC';
        end
        %Load GPFA and Decoder Params
        GPFAParams = load('D:\Animals\Earl\2020\08\20190821\analysis\mat_results\run005\gpfa_xDim10.mat')
        DecoderParams = BC(1).Decoder.Parameters;    
        

   
    HC = load('D:\Animals\Earl\2020\08\20190821\10_centerOut_HC\Earl20190821_10_centerOut_HC_SI_translated.mat');   
        %Get Data Struct
        HC = getDataStruct(HC.Data,'getMarker',true,'getSpikes',true); 
        numTrials = size(HC,2);
        for trial = 1:numTrials
           HC(trial).Task = 'HC';
        end
        
  Iso = load('D:\Animals\Earl\2020\08\20190821\08_isoForce\Earl20190821_08_isoForce_SI_translated.mat');
        %Get Data Struct
        Iso = getDataStruct(Iso.Data,'getMarker',true,'getSpikes',true); 
        numTrials = size(Iso,2);
        for trial = 1:numTrials
           Iso(trial).Task = 'Iso';
        end  
        

    %Combine data, keep only successful trials
    Data = [BC,HC,Iso];
    Data = Data([Data.trialStatus]==1);
    clearvars BC HC Iso
    

%% Label Conditions
    numTrials = size(Data,2);
    for trial = 1:numTrials
       Data(trial).Condition = [Data(trial).Task,'_T',num2str(Data(trial).targetData.targetID)]; 
    end

    
% %% Get minimum number of successful trials for any condition
%     conditionList = unique({Data.Condition});
%     minNumber = size(Data,2);
%     maxNumber = 0;
%     for condition = conditionList
%        tempData = Data(strcmpi({Data.Condition},condition))
%        trialStatus = [tempData.trialStatus];
%        numSuc = sum(trialStatus);
%        if numSuc < minNumber
%            minNumber = numSuc;
%        end
%        if numSuc > maxNumber
%            maxNumber = numSuc;
%        end
%     end
%     
%     
% %% Keep only minimum number
%     minSucData = [];
%     for condition = conditionList
%         tempData = Data(strcmpi({Data.Condition},condition));
%         minSucData = [minSucData,tempData(1:minNumber)];
%     end
%     Data = minSucData;
%     clearvars minSucData


end