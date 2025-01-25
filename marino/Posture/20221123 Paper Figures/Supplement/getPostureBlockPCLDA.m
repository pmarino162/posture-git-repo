function [posturePCLDA] = getPostureBlockPCLDA(trajStruct)

    % For each timecourse in trajstruct, get single point
    for i = 1:size(trajStruct,2)
        %Individual trials
        for j = 1:size(trajStruct(i).allSmoothFR,2)
           trajStruct(i).allSmoothFR(j).obs = mean(trajStruct(i).allSmoothFR(j).traj); 
        end
    end
        
    % Do PCA on observations 
    allObs = [];
    for i = 1:size(trajStruct,2)
        for j = 1:size(trajStruct(i).allSmoothFR,2)
            obs = trajStruct(i).allSmoothFR(j).obs;
            allObs = vertcat(allObs,obs);
        end
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allObs); 
    
    % Get posture list
    postureList = unique([trajStruct.posture]);
    
    %Form obsStruct
    obsStruct = struct('label',[],'numObs',[],'allObs',[]);
    structInd = 1; 
    for posture = postureList
        obsStruct(structInd).label = posture;
        allObs = [];
        tempTrajStruct = trajStruct([trajStruct.posture]==posture);
        for i = 1:size(tempTrajStruct,2)
            for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                obs = mean(tempTrajStruct(i).allSmoothFR(j).traj,1)*coeff(:,1:20);
                allObs = vertcat(allObs,obs);
            end
        end
        obsStruct(structInd).allObs = allObs;
        obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
        structInd = structInd + 1;
    end
    
    %Do LDA on Observations
        %Preallocate
        totalNumObs = sum([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = NaN(totalNumObs,numDims); labels =  NaN(totalNumObs,1); 

        %Fill
        k = 1;
        for i = 1:size(obsStruct,2)
            totalClassObs = size(obsStruct(i).allObs,1);
            obs(k:k+totalClassObs-1,:) = obsStruct(i).allObs;
            labels(k:k+totalClassObs-1,1) = ones(totalClassObs,1).*obsStruct(i).label;
            k = k+totalClassObs;
        end

        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
        
    %Combine PCA and LDA
    posturePCLDA = coeff(:,1:20)*LDAproj;
    
end

