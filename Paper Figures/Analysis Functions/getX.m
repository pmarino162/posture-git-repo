function [X] = getX(trajStruct,field,minNumTimestamps,postureList,numPostures,targetList,numTargets)    
    
        %Get 'X', which is the name used for tensor of de-meaned concatentated
        %condition-averages organized by posture and target. Field
        %arguments specifies whether to perform on z-scored FRs or PCA
        %scores.
        
        %Get numDims based on type of data
        switch field
            case 'avgZSmoothFR'
                numDims = size(trajStruct(1).avgZSmoothFR(1).traj,2);
            case 'avgPCA'
                numDims = size(trajStruct(1).avgPCA(1).traj,2);
            case 'avgSmoothFR'
                numDims = size(trajStruct(1).avgSmoothFR(1).traj,2);
        end
        
        %Get X
        X = NaN(minNumTimestamps,numTargets,numPostures,numDims);
        postureInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                switch field
                    case 'avgZSmoothFR'
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.traj;
                    case 'avgPCA'
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgPCA.traj;
                    case 'avgSmoothFR'
                        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
                end
                X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        
        %De-mean (for clarity, written to get the mean for each neuron
        %across conditions, then subtract)
        Xmean = mean(X,[1 2 3],'omitnan');
        X = X-Xmean;
end