clear; clc;

%% Load dataset
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    
%% Format spikes for neural traj  
    numTrials = size(Data,2);
    dat = struct('trialId',[],'spikes',[]);
    for trial = 1:numTrials
        dat(trial).trialId = trial;
        dat(trial).spikes = Data(trial).spikes(1:20,:);
    end 
    
%% Run neural traj
    result = neuralTraj('run000', dat);
    
%% Run causalInferencePrecomp
    
    precomp = causalInferencePrecomp(estParams, 22)