function [predictedLabels,pctCorrect] = NBClassify(observations,labels,numClasses,NBModel) 

% %Observations is NxD, where N is number of observations, and D is the
% %dimensionality of those observations
% 
% %Labels is Nx1
% %Author: Patrick Marino 3/8/2022
    
    %Predict labels 
    numObs = size(observations,1);
    numDims = size(observations,2);
    predictedLabels = zeros(1,numObs);
    for obs = 1:numObs
       obsData = observations(obs,:);
       classObjFcn = zeros(1,numClasses);
       for classInd = 1:numClasses
           mu = NBModel(classInd).mean;
           sigma = NBModel(classInd).var;
           for dim = 1:numDims
               classObjFcn(classInd) = classObjFcn(classInd) - (((obsData(dim)-mu(dim)).^2)/sigma(dim).^2);
           end
       end
       [~,predictedClassInd] = max(classObjFcn);
       predictedLabels(obs) = NBModel(predictedClassInd).label;
    end
    
    %Evaluate classification performance
    pctCorrect = 100*(sum(labels'==predictedLabels)./size(labels,1));
        
end