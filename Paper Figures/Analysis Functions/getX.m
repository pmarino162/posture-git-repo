function [X] = getX(trajStruct,minNumTimestamps,postureList,numPostures,targetList,numTargets,numChannels)    
    
        %Get 'X', which is the name used for tensor of concatentated
        %condition-averages organized by posture and target 
        X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
        postureInd = 1;
        for posture = postureList
            targetInd = 1;
            for target = targetList
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgZSmoothFR.traj;
                X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        
end