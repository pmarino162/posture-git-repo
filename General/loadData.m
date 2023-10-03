function [Data,zScoreParams] = loadData(dataset)

%Loads and returns a dataset

saveDir = getExperimentSaveDir(dataset);
load(fullfile(saveDir,[dataset,'.mat']))    

end