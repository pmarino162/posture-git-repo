function [kinStruct] = getKinStruct(Data,condFields,kinFields)

%% Get list of condition labels
    numTrials = size(Data,2);
    numCondFields = size(condFields,2);
    condLabels = zeros(numTrials,numCondFields);
    for trial = 1:numTrials
        for condField = 1:numCondFields
             condLabels(trial,condField) = getfield(Data(trial),condFields{condField}{2:end});
        end
    end
    condList = unique(condLabels,'rows');
    numCond = size(condList,1);
    
%% Create kinStruct, preallocate
    kinStruct = struct();
    for condField = 1:numCondFields
         condFieldName = condFields{condField}{1};
         kinStruct.(condFieldName) = [];
    end
    numKinFields = size(kinFields,2);
    for kinField = 1:numKinFields
         kinFieldName = kinFields{kinField};
         upperKinFieldName = [upper(kinFieldName(1)),kinFieldName(2:end)];
         kinStruct.(['all',upperKinFieldName]) = [];
    end
    kinStruct = repmat(kinStruct,1,numCond);
    
%% Fill kinStruct
    structInd = 1;
    for condInd = 1:numCond
        %Fill in condition information in kinStruct
        cond = condList(condInd,:);
        for condField = 1:numCondFields
            condFieldName = condFields{condField}{1};
            kinStruct(structInd).(condFieldName) = cond(condField);
        end
        %Get condition data
        tempData = Data(all(condLabels==cond,2));
        numCondTrials = size(tempData,2);        
        %Store all kin data
        for trial = 1:numCondTrials
            for kinField = 1:numKinFields
               kinFieldName = kinFields{kinField};
               upperKinFieldName = [upper(kinFieldName(1)),kinFieldName(2:end)];
               kindata = tempData(trial).kinData.(kinFieldName);
               if ~isempty(kindata)
                    kinStruct(structInd).(['all',upperKinFieldName])(trial) = kindata;
               else
                   kinStruct(structInd).(['all',upperKinFieldName])(trial) = NaN;
               end
            end
        end        
        %Update structInd
        structInd = structInd + 1;
    end

end