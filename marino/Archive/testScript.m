clear;
clc;
% clf;

%% Load Data
% temp1 = load('S:\Animals\Ford\2019\07\20190703\05_centerOut_03\Ford20190703_05_centerOut_03_SI_SW_translated.mat') 
% % 
% % %temp1 = load('S:\Animals\Ford\2019\06\20190624\03_centerOut_01\Ford20190624_03_centerOut_01_SI_SW_translated.mat');
%  
% % temp1 = load('D:\Animals\Earl\2019\04\20190429\06_centerOut_BC_memoryGuided_HC_01\Earl20190429_06_centerOut_BC_memoryGuided_HC_01_SI_translated.mat');
% rawData = [temp1.Data];
% clearvars temp1
animal = 'Dwight';
date = '20190806';
rawData = loadBCData(animal,date);
inclStates = {'BC Freeze','Touch Bar BC'};

%% Process Data 
Data = getDataStruct(rawData,'getWaveforms',true);
% clearvars rawData

%% Generate Success Rate Plot
plotSuccess(Data)

%% Plot Tuning Data
startState = 'Touch Bar BC';
windowBounds = [-100 200];
plotTuning(Data,startState,windowBounds);


%% Plot Drift
plotDrift(Data,inclStates,date);