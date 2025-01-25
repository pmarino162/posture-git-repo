clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis\reachVsBCI';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Construct resultStruct for multiple sessions
%     %Earl Within Day Reach vs BCI
%     inputStruct(1).dataset = 'E20200314'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'E20200314'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'E20200314'; inputStruct(3).task = 'Iso'; inputStruct(3).epoch = 'all';

    %Earl across days reach vs BCI
    inputStruct(1).dataset = 'E20200316'; inputStruct(1).task = 'BCI'; inputStruct(1).epoch = 'all';
    inputStruct(2).dataset = 'E20200317'; inputStruct(2).task = 'BCI'; inputStruct(2).epoch = 'all';
    inputStruct(3).dataset = 'E20200318'; inputStruct(3).task = 'BCI'; inputStruct(3).epoch = 'all';
    inputStruct(4).dataset = 'E20210706'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'all';
    inputStruct(5).dataset = 'E20210707'; inputStruct(5).task = 'Reach'; inputStruct(5).epoch = 'all';
    inputStruct(6).dataset = 'E20210708'; inputStruct(6).task = 'Reach'; inputStruct(6).epoch = 'all';
    
    %Nigel reach vs BCI
%     inputStruct(1).dataset = 'N20190222'; inputStruct(1).task = 'Reach'; inputStruct(1).epoch = 'all';
%     inputStruct(2).dataset = 'N20190226'; inputStruct(2).task = 'Reach'; inputStruct(2).epoch = 'all';
%     inputStruct(3).dataset = 'N20190227'; inputStruct(3).task = 'Reach'; inputStruct(3).epoch = 'all';
%     inputStruct(4).dataset = 'N20190228'; inputStruct(4).task = 'Reach'; inputStruct(4).epoch = 'all';
%     inputStruct(5).dataset = 'N20171215'; inputStruct(5).task = 'BCI'; inputStruct(5).epoch = 'all';
%     inputStruct(6).dataset = 'N20180221'; inputStruct(6).task = 'BCI'; inputStruct(6).epoch = 'all';

    %Rocky reach vs BCI



    resultStruct = [];
    for i = 1:numel(inputStruct)
        dataset = inputStruct(i).dataset;
        task = inputStruct(i).task;
        epoch = inputStruct(i).epoch;
        [tempResultStruct] = leaveOneConditionOutCompAnalysisReachVsBCI(dataset,task,epoch);
        resultStruct = vertcat(resultStruct,tempResultStruct);
    end

%% Get monkey info for plotting
    monkey = resultStruct(1).dataset(1);
    
 %% Create R2 Histogram 
figure ; hold on
    task = {resultStruct.task};
    histogram([resultStruct(strcmpi(task,'BCI')).R2],'BinWidth',.1);
    histogram([resultStruct(strcmpi(task,'Reach')).R2],'BinWidth',.1);
    xlim([0 1])
    xlabel('R^2')
    ylabel('Num Sessions')
    legend({'BCI','Reach'})  
    
%% Create mean dist histogram 
figure ; hold on
    task = {resultStruct.task};
    histogram([resultStruct(strcmpi(task,'BCI')).fullDist],'BinWidth',.1);
    histogram([resultStruct(strcmpi(task,'Reach')).fullDist],'BinWidth',.1);
    xlabel('Mean Dist (std)')
    ylabel('Num Sessions')
    legend({'BCI','Reach'})  
    [h,p] = ttest2([resultStruct(strcmpi(task,'BCI')).fullDist],...
        [resultStruct(strcmpi(task,'Reach')).fullDist])
    p = round(p,3)
    ax = gca;
    xlims = ax.XLim;
    ylims = ax.YLim;
    text(xlims(1)+.1,ylims(2)-.1,['p=',num2str(p)])
    if saveFig
        saveas(gcf,fullfile(saveDir,[monkey,'_','reachVsBCIMeanDist.svg']));
    end
    