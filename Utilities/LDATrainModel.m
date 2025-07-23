function [LDAModel,numClasses] = LDATrainModel(observations,labels)

% %Observations is NxD, where N is number of observations, and D is the
% %dimensionality of those observations
% 
% %Labels is Nx1
% %Author: Patrick Marino 2/14/2021

    %Define output 
    LDAModel = struct('label',[],'prior',[],'mean',[],'estCov',[]);
    
    %Get number of classes and observation dimensionality
    uniqueLabels = unique(labels)';
    numClasses = size(uniqueLabels,2);
    numDims = size(observations,2);
    
    %Estimate class covariances
    estCov = zeros(numDims);
    %Average across classes
    for class = uniqueLabels
       classData = observations(labels==class,:);
       classCov = cov(classData);
       estCov = estCov + classCov;
    end
    estCov = estCov./numClasses;
    
    %Get class priors (assume uniform distribution), means, and estimated
    %covariances
    classInd = 1;
    for class = uniqueLabels
       classData = observations(labels==class,:);
       classPrior = 1./numClasses; 
       classMean = mean(classData,1);
       LDAModel(classInd).label = class;
       LDAModel(classInd).prior = classPrior;
       LDAModel(classInd).mean = classMean;
       LDAModel(classInd).estCov = estCov;
       classInd = classInd + 1;
    end
    
    %Estimate class means and covariances
%     classMeans = zeros(numClasses,numDims);
%  classMeans(class,:) = mean(classData,1);
%     classPriors = zeros(1,numClasses);
end