function [Data,zScoreParams] = preprocessAndSaveData20220419(dataset)

%Loads raw dataset, preprocesses it, gets Data struct, and saves preprocessed struct. Input is a string
%that specifies the monkey and date. Each dataset has a specific function
%for implementing dataset-specific instructions.

%% Get Data Struct
    switch dataset
        %% Earl
        case {'E20200316','E20200317','E20200318','E20200319'} %Earl 5-posture BCI CO: E30:I30
            %[Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
            [Data,~,~,~] = loadEarlData_5Posture_BCI_CO(dataset);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case {'E20210901'} %Earl Shoulder + Elbow Rotation BCI
            [Data] = loadEarlData20210901_20211210;
            [Data] = binSpikes1msAddTime(Data);
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'} %Earl 7-posture DCO
            %[Data] = loadEarlData20210706_20211210;
            [Data] = loadEarlData_7postureReaching(dataset);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case {'E20200311','E20200312','E20200313','E20200314'} %Earl Multiple tasks paradigm
            %[Data,~,~,~,~,~,~] = loadEarlData20200314_20211210;
            [Data] = loadEarlData_3Posture_MultipleTask(dataset);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
            %Remove reaches that didn't begin with hand on forcebar
            [Data] = removeBadMultipleTasksReachingTrials(Data);            
        case {'E20200116','E20200117','E20200120'} %Earl 5-posture isometric force
            %[Data] = loadEarlData20200116_20211210();
            [Data,~] = loadEarlData_5Posture_IsoForce(dataset);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case 'E20211007' %Earl 8-target isoforce; intrinsic rotation
            [Data] = loadEarlData20211007_20211210();
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case 'E20190830' %Earl 5 posture CO for drift control
            [Data] = loadEarlData_20190830;
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case 'E20190729'%'E20190728' %Earl 'alternate decoder' experiments
            %[Data] = loadEarlData20190728;
            [Data] = loadEarlData20190729;
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case 'E20210823' %Earl multiple decoders in 1 posture control
            [Data] = loadEarlData20210823_20211210;
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
            
            
        %% Nigel
        case {'N20171215','N20180221'} %Nigel BCI: N00, A90, and I45
            if strcmpi(dataset,'N20171215')
                [Data] = loadNigelData20171215_20220419;
            elseif strcmpi(dataset,'N20180221')
                [Data] = loadNigelData20180221_20220419;
            end
            %Standardize spikes struct across trials to include all sorts/units
            [Data] = standardizeSpikesStruct(Data);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
            
        case {'N20190222','N20190226','N20190227'... %Nigel Reaching: Dissociation Experiment
                'N20190228','N20190301','N20190305','N20190306','N20190307'}
            %[Data] = loadNigelData20190226_20220419;
            [Data,~] = loadChaseLabData_Dissociation(dataset);
            %Keep only lower-right visual data
            conditionData = [Data.conditionData];
            Data = Data([conditionData.visualID]==2);
            %Replicate time field in marker field so that format matches Earl's
            for trial = 1:numel(Data)
                Data(trial).marker.time = Data(trial).time;
            end
            
        %% Rocky 
        case 'R20201020' %Rocky BCI: N00 and A90
            %Convert from intermediate to formatted (if necessary)
                %int2formatted_Marino(dataset);
            %Load/preprocess formatted data
            [Data,postureIDs] = loadRockyData20201020_20220419;
            %Standardize spikes struct across trials to include all sorts/units
            [Data] = standardizeSpikesStruct(Data);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case 'R20201021' %Rocky BCI: N00 and I45
            %Convert from intermediate to formatted (if necessary)
                %int2formatted_Marino(dataset);
            %Load/preprocess formatted data
            [Data,postureIDs] = loadRockyData20201021_20220419;
            %Standardize spikes struct across trials to include all sorts/units
            [Data] = standardizeSpikesStruct(Data);
            %Bin spikes at 1ms, add time field corresponding to bins
            [Data] = binSpikes1msAddTime(Data);
        case {'R20200221','R20200222'} %Rocky Reaching: Dissociation Experiment 
            [Data,~] = loadChaseLabData_Dissociation(dataset);
            %Keep only lower-right visual data
            conditionData = [Data.conditionData];
            Data = Data([conditionData.visualID]==2);
            %Replicate time field in marker field so that format matches Earl's
            for trial = 1:numel(Data)
                Data(trial).marker.time = Data(trial).time;
            end
            %Clean data (remove a few select bad trials)
            if strcmp(dataset,'R20200221')
                [Data] = cleanDataR20200221(Data);
            end
    end

%% Run Common Preprocessing Scripts
    %Keep only successful trials
    switch dataset
        case 'E20190728'
        otherwise
            Data = Data([Data.trialStatus]==1);
    end
    %Remove "no-spike" trials
        %Data = removeNoSpikeTrials(Data);
        
    %Filter Marker
    
    
    %Get Kinematic Data
        [Data] = getKinData20220419(Data);
        
    %Apply trial exclusion criteria to main datasets
    switch dataset
        case 'E20190729'
        otherwise
            [Data] = applyExclusionCriteria(Data);
    end
    %Get z-score or softnorm params. Use sortMeans from 25ms bin width to
    %remove Low FR Sorts
        structInd = 1;
        for binWidth = [25,50]
           kernelStdDev = binWidth; 
           [sortMean,sortStd,sortRange] = getZScoreParams(Data,binWidth,kernelStdDev);
           if binWidth == 25
                threshold = 3;
                rmSorts = sortMean < threshold;
                numTrials = size(Data,2);
                for trial = 1:numTrials
                    Data(trial).spikes(:,rmSorts) = [];
                end
                sortMean = sortMean(~rmSorts);
                sortStd = sortStd(~rmSorts);
                sortRange = sortRange(~rmSorts); 
           end
           zScoreParams(structInd).binWidth = binWidth;
           zScoreParams(structInd).kernelStdDev = kernelStdDev;
           zScoreParams(structInd).sortMean = sortMean;
           zScoreParams(structInd).sortStd = sortStd;
           zScoreParams(structInd).sortRange = sortRange; 
           structInd = structInd + 1;
        end
    
    
    %Remove coincident channels    
    

     
    
%% Save Preprocessed Dataset    
    saveDir = getExperimentSaveDir20220419(dataset);
    save(fullfile(saveDir,[dataset,'.mat']),'Data','zScoreParams','-v7.3');

end