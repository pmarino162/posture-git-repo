clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Figure 3\Projections';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
        load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\purpleAndGreen.mat')
    tcmap = purpleAndGreen;

%% Run loop for each dataset      
    %Load data
    dataset = 'R20200221';
    [Data,zScoreParams] = loadData(dataset);

   
%% Get rxnTimes

rxnTime = [Data.kinData];
rxnTime = [rxnTime.rxnTime];
figure
histogram(rxnTime);