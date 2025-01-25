clear; clc; clf; close all

%% Setup save directory
    saveDir = 'D:\Posture\Joel Data Transfer';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Load Earl reaching; eliminate fields; save
    %'E20210706','E20210707','E20210708','E20210709','E20210710'
    dataset = 'E20210710';
    [Data,zScoreParams] = loadData(dataset);
    Data = rmfield(Data,{'force','Decoder'});
    save(fullfile(saveDir,['DelayedCenterOut_',dataset,'.mat']),'Data','-v7.3');

%% Load Earl BCI; eliminate fields; save
    %'E20200316','E20200317','E20200318','E20200319'
    dataset = 'E20200319'; 
    [Data,zScoreParams] = loadData(dataset);
    Data = rmfield(Data,{'force','kinData'});
    save(fullfile(saveDir,['BCI_',dataset,'.mat']),'Data','-v7.3');