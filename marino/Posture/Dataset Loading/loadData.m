function [Data,zScoreParams] = loadData(dataset)

%Loads and returns a dataset

saveDir = getExperimentSaveDir20220419(dataset);
load(fullfile(saveDir,[dataset,'.mat']))    

end