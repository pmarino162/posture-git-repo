function [allTraj] = collectAllAvgTraj(trajStruct)
        
        %Get 'allTraj', which is the name used for matrix of concatentated
        %condition-averages (mean-centered, truncated to minimum length)
        
        [~,~,~,~,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct);

        allTraj = NaN(numConditions*minNumTimestamps,numChannels);
        j = 1;
        for i = 1:numConditions
           allTraj(j:j+minNumTimestamps-1,:) = trajStruct(i).avgZSmoothFR.traj(1:minNumTimestamps,:);
           j = j + minNumTimestamps;
        end
        
        %Mean center
        allTraj = allTraj - mean(allTraj,1);
        
end