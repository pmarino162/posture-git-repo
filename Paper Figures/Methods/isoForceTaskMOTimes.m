clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Figure 3\Projections';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Run loop for each dataset      
    %Load data
    dataset = 'E20200116';
    [Data,zScoreParams] = loadData(dataset);
   
%% Get rxnTimes
    kinData = [Data.kinData];
    rxnTime = [kinData.rxnTime];
    figure
    histogram(rxnTime);