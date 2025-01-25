clear; clc; clf; close all;
    
%% Setup save dir   
    dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210930 - first model\';

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = rainbow2d;
    
%% Load Data 
    binWidth = 25;
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    exclCh =  [44 87 88 77 78 71 67 69 118];
    getSorts = false;
    Data = getDataStruct(Data,'getSpikes',true,'getSorts',getSorts,'exclCh',exclCh,'binWidth',binWidth,'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
    allTrialPostures = [Data.conditionData];
    allTrialPostures = [allTrialPostures.postureID];
%     keepPostures = [1,2,4,5];
    keepPostures = [1:7];
    Data = Data(ismember(allTrialPostures,keepPostures));
    
%% Plot target posture map
    conditionData = [Data.conditionData];
    targetData = [Data.targetData];
    figure
    hold on
    for posture = [5]
        tempData = Data([conditionData.postureID]==posture & [targetData.targetID]==1);
        workspaceCenter = tempData(1).targetData.workspaceCenter;
%         text(workspaceCenter(1),workspaceCenter(2),num2str(posture),'Color',tcmap(posture,:),'FontSize',14) 
        if ismember(posture,[5])
            for target = 1:8
                tempData = Data([conditionData.postureID]==posture & [targetData.targetID]==target);
                if ~isempty(tempData)
                    targetLoc = tempData(1).targetData.targetLoc;
                    targetLoc = targetLoc + workspaceCenter;
                    plot(targetLoc(1),targetLoc(2),'.','MarkerSize',20,'Color',tcmap(target,:));
                end
            end
        end
    end
%     set(gca,'visible','off')
    axis off
    ax = gca;
    ax.XLim = [-225 200];
    ax.YLim = [-625 -200];
    xlabel('x (mm)')
    ylabel('y (mm)')
