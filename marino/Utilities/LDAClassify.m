function [predictedLabels,pctCorrect] = LDAClassify(observations,labels,numClasses,LDAModel) 

% %Observations is NxD, where N is number of observations, and D is the
% %dimensionality of those observations
% 
% %Labels is Nx1
% %Author: Patrick Marino 2/14/2021
    
    %Predict labels 
    numObs = size(observations,1);
    predictedLabels = zeros(1,numObs);
    for obs = 1:numObs
       obsData = observations(obs,:);
       classObjFcn = zeros(1,numClasses);
       for classInd = 1:numClasses
           prior = LDAModel(classInd).prior;
           mu = LDAModel(classInd).mean;
           sigma = LDAModel(classInd).estCov;
%           prior = classPriors(class);
%           mu = classMeans(class,:);
%           sigma = estCov;
          classObjFcn(classInd) = log(prior) - 0.5*(obsData-mu)*inv(sigma)*(obsData-mu)';
       end
       [~,predictedClassInd] = max(classObjFcn);
       predictedLabels(obs) = LDAModel(predictedClassInd).label;
%        [~,predictedLabels(obs)] = max(classObjFcn);
    end
    
    %Evaluate classification performance
    pctCorrect = 100*(sum(labels'==predictedLabels)./size(labels,1));
        
end