function [NBModel,numClasses] = NBTrainModel(observations,labels)

% %Observations is NxD, where N is number of observations, and D is the
% %dimensionality of those observations
% 
% %Labels is Nx1
% %Author: Patrick Marino 3/8/2022

    %Define output 
    NBModel = struct('label',[],'mean',[],'var',[]);
    
    %Get number of classes and list of labels
    uniqueLabels = unique(labels)';
    numClasses = size(uniqueLabels,2);
    
    %Estimate class means and variances for each neuron (assume same variance across
    %classes
    classInd = 1;
    for class = uniqueLabels
       classData = observations(labels==class,:);
       classMean = mean(classData,1);
       NBModel(classInd).label = class;
       NBModel(classInd).mean = classMean;
       NBModel(classInd).var = var(observations);
       classInd = classInd + 1;
    end

end