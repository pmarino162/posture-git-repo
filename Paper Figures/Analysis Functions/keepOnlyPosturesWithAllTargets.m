function [trajStruct] = keepOnlyPosturesWithAllTargets(trajStruct,postureList,targetList)

    keepPosture = [];
    for posture = postureList
        tempTrajStruct = trajStruct([trajStruct.posture]==posture);
        postureTargetList = [tempTrajStruct.target];
        if isequal(postureTargetList,targetList)
            keepPosture = [posture,keepPosture];
        end
    end
    trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));

end