function [trajStruct,PCs] = projectToTopPCs(allTraj,trajStruct,numPCsToKeep)

        %Reduces dimensionality of 'allTraj', keeping only top PCs
        %allTraj is already mean-centered, so no need to mean-center again
        %before projection

        [PCs,~,~,~,~,~] = pca(allTraj);
        
        for i = 1:size(trajStruct,2)
           trajStruct(i).avgPCA.traj = trajStruct(i).avgZSmoothFR.traj*PCs(:,1:numPCsToKeep); 
           for j = 1:size(trajStruct(i).allZSmoothFR,2)
                trajStruct(i).allPCA(j).traj = trajStruct(i).allZSmoothFR(j).traj*PCs(:,1:numPCsToKeep);
                trajStruct(i).allPCA(j).timestamps = trajStruct(i).allZSmoothFR(j).timestamps;
           end
        end
        
        
end