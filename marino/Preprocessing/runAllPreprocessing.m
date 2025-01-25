%% Reaching
%Earl
    for dataset = {'E20210706','E20210707','E20210708','E20210709','E20210710'}
        [Data,zScoreParams] = preprocessAndSaveData20220419(dataset{1,1});
    end

%Rocky
    for dataset = {'R20200221','R20200222'}
        [Data,zScoreParams] = preprocessAndSaveData20220419(dataset{1,1});
    end
%     
    %Nigel
    for dataset = {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190306','N20190307'}
        [Data,zScoreParams] = preprocessAndSaveData20220419(dataset{1,1});
    end