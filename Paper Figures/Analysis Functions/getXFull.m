function [Xfull] = getXFull(trajStruct,field,minNumTimestamps,postureList,numPostures,targetList,numTargets)    
    
        %Get 'Xfull', which is the name used for tensor of de-meaned concatentated
        %trials organized by posture and target. Field
        %arguments specifies whether to perform on z-scored FRs or PCA
        %scores.
        
        %Get numDims based on type of data
        switch field
            case 'allZSmoothFR'
                numDims = size(trajStruct(1).allZSmoothFR(1).traj,2);
            case 'allPCA'
                numDims = size(trajStruct(1).allPCA(1).traj,2);
            case 'allSmoothFR'
                numDims = size(trajStruct(1).allSmoothFR(1).traj,2);
        end
        
        %Get max num cond trials; preallocate X
        [maxNumCondTrials] = getMaxNumCondTrials(trajStruct);
        Xfull = NaN(minNumTimestamps,numTargets,numPostures,numDims,maxNumCondTrials);
        
        %Get Xfull        
        postureInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                %Get numCondTrials
                switch field
                    case 'allZSmoothFR'
                        numCondTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR,2);
                    case 'allPCA'
                        numCondTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allPCA,2);
                    case 'allSmoothFR'
                        numCondTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR,2);
                end
                %For each cond trial, store traj
                for trial = 1:numCondTrials
                    switch field
                        case 'allZSmoothFR'
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allZSmoothFR(trial).traj;
                        case 'allPCA'
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allPCA(trial).traj;
                        case 'allSmoothFR'
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR(trial).traj;
                    end
                    Xfull(:,targetInd,postureInd,:,trial) = traj(1:minNumTimestamps,:);     
                end
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        
        %De-mean (for clarity, written to get the mean for each neuron
        %across conditions, then subtract)
        Xmean = mean(Xfull,[1 2 3 5],'omitnan');
        Xfull = Xfull-Xmean;
end