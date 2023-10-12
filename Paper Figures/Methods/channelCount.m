clear; clc; clf; close all

%% Main loop
    resultStruct = struct('monkey',[],'allChCount',[],'uChCount',[],'stdChCount',[]);
    structInd = 1;
    for monkey = {'E1'}%,'E2','N','R'}

        %Get datasetList
        switch monkey{1}
            case 'E1'
                datasetList = {'E20200316','E20200317','E20200318','E20200319',...
                    'E20200116','E20200117','E20200120'};
            case 'E2'
                datasetList = {'E20210706','E20210707','E20210708','E20210709','E20210710'};
            case 'N'
                datasetList = {'N20171215','N20180221',...
                    'N20190222','N20190226','N20190227','N20190228','N20190301','N20190307'};
            case 'R'
                datasetList = {'R20201020','R20201021',...
                    'R20200221','R20200222'};
        end

        %Get allChCount
        numDatasets = length(datasetList);
        allChCount = zeros(1,numDatasets);
        for datasetInd = 1:numDatasets
            dataset = datasetList{datasetInd};
            [Data,zScoreParams] = loadData(dataset);
            allChCount(datasetInd) = size(Data(1).spikes,2);
        end

        %Store result
        resultStruct(structInd).monkey = monkey;
        resultStruct(structInd).allChCount = allChCount;
        resultStruct(structInd).uChCount = mean(allChCount);
        resultStruct(structInd).stdChCount = std(allChCount);
        structInd = structInd + 1;
    end