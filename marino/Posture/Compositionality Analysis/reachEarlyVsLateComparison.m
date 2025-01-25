clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis\reachEarlyVsLate';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Construct resultStruct for multiple sessions
%     %Rocky reaching early vs. late
%     inputStruct(1).dataset = 'R20200221'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'early';
%     inputStruct(2).dataset = 'R20200221'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'late';
%     inputStruct(3).dataset = 'R20200222'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'early';
%     inputStruct(4).dataset = 'R20200222'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'late';

    %Earl reaching early vs. late
    inputStruct(1).dataset = 'E20210706'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'early';
    inputStruct(2).dataset = 'E20210706'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'late';
    inputStruct(3).dataset = 'E20210707'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'early';
    inputStruct(4).dataset = 'E20210707'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'late';
    inputStruct(5).dataset = 'E20210708'; inputStruct(5).task = 'Reach'; inputStruct(5).epoch = 'early';
    inputStruct(6).dataset = 'E20210708'; inputStruct(6).task = 'Reach'; inputStruct(6).epoch = 'late';
    
%     %Nigel reaching early vs. late
%     inputStruct(1).dataset = 'N20190222'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'early';
%     inputStruct(2).dataset = 'N20190222'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'late';
%     inputStruct(3).dataset = 'N20190226'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'early';
%     inputStruct(4).dataset = 'N20190226'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'late';
%     inputStruct(5).dataset = 'N20190227'; inputStruct(5).task = 'Reach'; inputStruct(5).epoch = 'early';
%     inputStruct(6).dataset = 'N20190227'; inputStruct(6).task = 'Reach'; inputStruct(6).epoch = 'late';
%     inputStruct(7).dataset = 'N20190228'; inputStruct(7).task = 'Reach'; inputStruct(7).epoch = 'early';
%     inputStruct(8).dataset = 'N20190228'; inputStruct(8).task = 'Reach'; inputStruct(8).epoch = 'late';
    
    resultStruct = [];
    for i = 1:numel(inputStruct)
        dataset = inputStruct(i).dataset;
        task = inputStruct(i).task;
        epoch = inputStruct(i).epoch;
        [tempResultStruct] = leaveOneConditionOutCompAnalysis(dataset,task,epoch);
        resultStruct = vertcat(resultStruct,tempResultStruct);
    end
    
%% Get monkey info for plotting
    monkey = resultStruct(1).dataset(1);

 %% Create R2 Histogram 
figure ; hold on
    epoch = {resultStruct.epoch};
    histogram([resultStruct(strcmpi(epoch,'early')).R2],'BinWidth',.1);
    histogram([resultStruct(strcmpi(epoch,'late')).R2],'BinWidth',.1);
    xlim([0 1])
    xlabel('R^2')
    ylabel('Num Sessions')
    legend({'Early','Late'})  
    
%% Create mean dist histogram 
figure ; hold on
    epoch = {resultStruct.epoch};
    histogram([resultStruct(strcmpi(epoch,'early')).fullDist],'BinWidth',.1);
    histogram([resultStruct(strcmpi(epoch,'late')).fullDist],'BinWidth',.1);
    xlabel('Mean Dist (std)')
    ylabel('Num Sessions')
    legend({'Early','Late'})  
    [h,p] = ttest([resultStruct(strcmpi(epoch,'early')).fullDist],...
        [resultStruct(strcmpi(epoch,'late')).fullDist])
    p = round(p,3)
    ax = gca;
    xlims = ax.XLim;
    ylims = ax.YLim;
    text(xlims(1)+.1,ylims(2)-.1,['p=',num2str(p)])
    if saveFig
        saveas(gcf,fullfile(saveDir,[monkey,'_','EarlyVsLateMeanDist.svg']));
    end
