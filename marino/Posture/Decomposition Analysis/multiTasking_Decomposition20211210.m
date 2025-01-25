clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Load Data 
    [Data,I30GPFAParams,I30DecoderParams,...
     N00GPFAParams,N00DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200314_20211210;
 
%% Get trajStruct 
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    
    binWidth = 25;

    %Reach
    trialInclStates(1).trialName = {'GridReaching'};
    
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Success with Reward','first',0}};
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
    %BC
    
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'state','Target Hold','first',100},{'state','Success with Reward','first',0}};
    holdTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false); 