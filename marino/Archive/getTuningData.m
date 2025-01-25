function [tuningData,targetAngles] = getTuningData(Data,dataType)

    %% Get number of trials, targets, and neuralDims (channels, sorts, or factors)
        %Trials
        numTrials = size(Data,2);
        %Targets
        targetData = [Data.targetData];
        targetIDList = unique([targetData.target1ID]);
        numTargets = size(targetIDList,2);
        %neuralDims
        switch dataType
            case 'channels'
                traj = Data(1).Traj.spikesTraj;
            case 'sorts'
            case 'GPFATraj'
                traj = Data(1).Traj.GPFATraj;
        end
        numNeuralDims = size(traj,2);
    
    %% Get list of target angles
        targetAngles = zeros(1,numTargets);
        for targetID = targetIDList
            targetDataInd = min(find([targetData.target1ID]==targetID));
            centerLoc = targetData(targetDataInd).centerLoc;
            targetLoc =  targetData(targetDataInd).target1Loc;
            targetVec = (targetLoc - centerLoc)';
            targetVec(3,1) = 0;
            targetAng = atan2d(norm(cross([1;0;0],targetVec)),dot([1;0;0],targetVec));
            if targetVec(2,1) < 0
                targetAng = 360-targetAng;
            end
            targetAngles(targetID) = targetAng;
        end
        
    %% Preallocate tuningData
        allData = NaN(numTrials,numTargets);
        means = NaN(1,numTargets);
        SD = NaN(1,numTargets);
        b0 = NaN(1,numTargets);
        p = NaN(1,1);
        PD = NaN(1,1);
        MD = NaN(1,1);
        tuningData = struct('allData',allData,'means',means,'SD',SD,'b0',b0,'p',p,'PD',PD,'MD',MD);
        tuningData = repmat(tuningData,1,numNeuralDims);
        
    %% Collect relevant data from each trial.  Store in 'allData' field of tuningData.      
        for trial = 1:numTrials
            %Get target data
            targetID = Data(trial).targetData.target1ID;
            %Get neural data
            switch dataType
                case 'channels'
                    traj = Data(trial).Traj.spikesTraj;
                    binTimes = Data(trial).Traj.spikesTrajTimestamps;
                    if ~isempty(traj)
                        timeLength = binTimes(end)-binTimes(1);
                        chSum = sum(traj,1);
                        dataToStore = chSum./(timeLength/1000);
                    end
                case 'sorts'
                case 'GPFATraj'
                    traj = Data(trial).Traj.GPFATraj;
                    if ~isempty(traj)
                        dataToStore = mean(traj,1);
                    end
            end
            %Store neural data
            if ~isempty(traj)
                for dim=1:numNeuralDims
                    tuningData(dim).allData(trial,targetID) = dataToStore(dim);
                end
            end
        end
    
    %% Get minimum number of observations for any target
        numObs = zeros(1,numTargets);
        for target = 1:numTargets
            numObs(target) =  sum(double(~isnan(tuningData(1).allData(:,target))));
        end
        minNumObs = min(numObs);
        
    %% Fit tuning curves. Fit Means and SD's. Store in appropriate fields of tuningData
        for dim=1:numNeuralDims
            y = []; x = [];
            dimAllData = tuningData(dim).allData;
            for target = targetIDList
                %Get relevant data for this dim for this target
                dimTargetData = dimAllData(:,target);
                nanMask = isnan(dimTargetData);
                dimTargetData(nanMask) = [];
                dimTargetData = dimTargetData(1:minNumObs,1);
                %Get mean and SD
                tuningData(dim).means(1,target) = mean(dimTargetData);
                tuningData(dim).SD(1,target) = std(dimTargetData);
                %Add neural data and target data to regression matrices
                targData = [ones(minNumObs,1),ones(minNumObs,1)*sind(targetAngles(target)),ones(minNumObs,1)*cosd(targetAngles(target))];
                y = vertcat(y,dimTargetData);
                x = vertcat(x,targData);
            end
    %         y = tuningData(1).channels.Means(channel,:)';
                [B,bint,r,rint,stats] = regress(y,x);
    %             B = x\y; 
                b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
                tuningData(dim).MD = sqrt(b1.^2 + b2.^2);
                tuningData(dim).PD= atan2d(b1,b2);
    %         if atan2d(b1,b2) < 0
    %             tuningData.channels(channel).PD= 360+atan2d(b1,b2);
    %         end
            tuningData(dim).b0 = b0;
            tuningData(dim).p = p;
        end
        
    
end