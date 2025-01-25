function [actualTime] = getActualTime(analogChannelNames,analogData)
    %Get Time Ind
    memoryHolderInd = find(strcmpi(analogChannelNames,'Memory Holder for State ID'));
    analogChannelNames(memoryHolderInd) = [];
    timeInd = find(strcmpi(analogChannelNames,'ActualTime'));
    
    %Get Actual Time
    actualTime = double(analogData(:,timeInd))./2;
end