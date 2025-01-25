function [Data] = cleanData20210701(Data)

%Exclude Postures
    rmTrials = [];
    numTrials = size(Data,2);
    for trial = 1:numTrials
       WC = Data(trial).targetData.workspaceCenter(1:2); 
       if isequal(WC,[-22.5,-432.5]) || isequal(WC,[70,-455]) || isequal(WC,[20,-475]) || isequal(WC,[-65,-390]) || isequal(WC,[-130,-390]) 
           
       else
           rmTrials = [rmTrials,trial];
       end
    end
    Data(rmTrials) = [];
    
%Delete bad trials from -65,-390 posture
    targetData = [Data.targetData];
    WC = cell2mat({targetData.workspaceCenter}');
    allPostureTrials = find(sum(ismember(WC,[-65,-390,680]),2)==3);
    Data(allPostureTrials(1:38)) = [];
    
% %Combine Sorts & Keep only sorts >=2
%     combineSorts = [9 2 3; 15 2 1; 23 2 3; 26 2 1; 4 2 1; 7 2 3; 13 2 1;...
%         29 2 1; 42 2 1; 46 2 1; 54 2 3; 56 2 3];
% 
%     specialChannels = combineSorts(:,1)';
%     
%     numTrials = size(Data,2);
%     for trial = 1:numTrials
%         sortList = Data(trial).spikes.sortList;
%         finalSortList = sortList;
%         %Create final sort list 
%         for channel = specialChannels
%             rmSort = combineSorts(specialChannels==channel,3);
%             rmInd = find([finalSortList.channel]==channel & [finalSortList.sort]==rmSort);
%             finalSortList(rmInd) = [];
%         end
%         finalSortList([finalSortList.sort]<2) = [];
%         
%         %Create final allSortSpikeBins
%         allSortSpikeBins = Data(trial).spikes.allSortSpikeBins;
%         finalAllSortSpikeBins = zeros(size(allSortSpikeBins,1),size(finalSortList,2));
% 
%         finalSpikeBinsInd = 1;
%         for channel = 1:128
%            if ismember(channel,specialChannels)
%                 chSortList = sortList([sortList.channel]==channel & [sortList.sort]>=2);
%                 for i = 1:size(chSortList,2)
%                    sort = chSortList(i).sort;
%                    if sort == combineSorts([combineSorts(:,1)==channel],2)
%                       sortInd1 = find([sortList.channel]==channel & [sortList.sort]==sort);
%                       sortInd2 = find([sortList.channel]==channel & [sortList.sort]==combineSorts([combineSorts(:,1)==channel],3));
%                       finalAllSortSpikeBins(:,finalSpikeBinsInd) = allSortSpikeBins(:,sortInd1) +  allSortSpikeBins(:,sortInd2);
%                       finalSpikeBinsInd = finalSpikeBinsInd + 1;
%                    elseif sort == combineSorts([combineSorts(:,1)==channel],3)
%                        %Skip it
%                    else
%                       sortInd = find([sortList.channel]==channel & [sortList.sort]==sort);
%                       finalAllSortSpikeBins(:,finalSpikeBinsInd) = allSortSpikeBins(:,sortInd);
%                       finalSpikeBinsInd = finalSpikeBinsInd + 1;
%                    end
%                 end
%                 
%                 validSortInds = find([sortList.channel]==channel & [sortList.sort] >= 2);
%                 for i = validSortInds
%                     finalAllSortSpikeBins(:,finalSpikeBinsInd) = allSortSpikeBins(:,i);
%                     finalSpikeBinsInd = finalSpikeBinsInd + 1;
%                 end
%            else
%                 validSortInds = find([sortList.channel]==channel & [sortList.sort] >= 2);
%                 for i = validSortInds
%                     finalAllSortSpikeBins(:,finalSpikeBinsInd) = allSortSpikeBins(:,i);
%                     finalSpikeBinsInd = finalSpikeBinsInd + 1;
%                 end
%            end
%         end
%         
%         %Write Data
%         Data(trial).spikes.sortList = finalSortList;
%         Data(trial).spikes.allSortSpikeBins = finalAllSortSpikeBins;
%     end
end