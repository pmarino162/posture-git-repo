function [Data,postureIDs] = labelPostures20210706(Data)
    %% Label Postures
    gridReachData = Data(cellfun(@(x) strcmpi(x,'GridReaching'),{Data.trialName}));
    targetData = [gridReachData.targetData];
    workspaceCenter = cell2mat({targetData.workspaceCenter}');
    uniWorkspaceCenter = unique(workspaceCenter,'rows');
    sortedYs = sort(unique(uniWorkspaceCenter(:,2)),'descend')';
    sortedXs = sort(unique(uniWorkspaceCenter(:,1)))';
    numCol = size(sortedXs,2);
    numTrials = size(gridReachData,2);
    for trial = 1:numTrials
       workspaceCenter = gridReachData(trial).targetData.workspaceCenter(1:2);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
        if row == 1 && col == 1
            postureID = 1;
        elseif row == 1 && col == 2
            postureID = 2;
        elseif row == 1 && col == 4
            postureID = 3;
        elseif row == 1 && col == 5
            postureID = 4;
        elseif row == 2 && col == 3
            postureID = 5;
        elseif row == 3 && col == 2
            postureID = 6;
        elseif row == 3 && col == 4
            postureID = 7;
        end
        gridReachData(trial).conditionData.postureID = postureID;
       %gridReachData(trial).conditionData.postureID = (row-1)*numCol + col;
    end
    %Create postureIDs list
    postureIDs = struct('postureID',[],'workspaceCenter',[]);
    for i = 1:size(uniWorkspaceCenter,1)
       workspaceCenter = uniWorkspaceCenter(i,:);
       row = find(sortedYs==workspaceCenter(2));
       col = find(sortedXs==workspaceCenter(1));
       
         if row == 1 && col == 1
            postureID = 1;
        elseif row == 1 && col == 2
            postureID = 2;
        elseif row == 1 && col == 4
            postureID = 3;
        elseif row == 1 && col == 5
            postureID = 4;
        elseif row == 2 && col == 3
            postureID = 5;
        elseif row == 3 && col == 2
            postureID = 6;
        elseif row == 3 && col == 4
            postureID = 7;
        end
       
%        postureID = (row-1)*numCol + col;
       postureIDs(i).postureID = postureID;
       postureIDs(i).workspaceCenter = workspaceCenter;      
    end
    [~,sortInd] = sort([postureIDs.postureID]);
    postureIDs = postureIDs(sortInd);
    
    %% Keep only gridReachingData; 
    Data = gridReachData;
    clearvars gridReachData 

end