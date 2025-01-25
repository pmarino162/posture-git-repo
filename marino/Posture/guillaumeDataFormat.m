clear; clc; clf; close all

%% Setup save directory
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Guillaume Data Transfer\';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Load reaching; eliminate fields; save
    dataset = 'E20210706';
    [Data,zScoreParams] = loadData(dataset);
    Data = rmfield(Data,{'force','Decoder'});
    save(fullfile(saveDir,['DelayedCenterOut_',dataset,'.mat']),'Data','-v7.3');

%% Load BCI; eliminate fields; save
    dataset = 'E20200317';
    [Data,zScoreParams] = loadData(dataset);
    Data = rmfield(Data,{'force','kinData'});
    save(fullfile(saveDir,['BCI_',dataset,'.mat']),'Data','-v7.3');