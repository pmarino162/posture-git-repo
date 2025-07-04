function [Data] = subtractReachWorkspaceCenters(Data,dataset)

    numTrials = size(Data,2);
    
    for trial = 1:numTrials

        if ismember(dataset,{'E20210706','E20210707','E20210708'})
            workspaceCenter = Data(trial).targetData.workspaceCenter;
            Data(trial).marker.position = Data(trial).marker.position(:,1:2) - workspaceCenter(1,1:2);
        end

        if ismember(dataset,{'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307','R20200221', 'R20200222'})
            workspaceCenter = Data(trial).targetData.workspaceCenter;
            Data(trial).cursor.position = Data(trial).cursor.position - workspaceCenter;
        end
        
        
        
        

    end


end