clear; clc; close all;
    
%% Set up figure
    f=figure; hold on;

%% Define color map
    color = [249, 141, 249]./255;
    %color = [253, 187, 253]./255;
    numRows = 50;
    col1 = linspace(1, color(1), numRows);
    col2 = linspace(1, color(2), numRows);
    col3 = linspace(1, color(3), numRows);
    cmap = [col1',col2',col3'];
    
%% Shading 
    numRows = 200;
    numCols = 200;
    C = zeros(numRows,numCols);
    for col = 2:numCols
        numRowsBelowDiag = col-1;
        shadeIncr = 1/numRowsBelowDiag;
        for row = (numRows-(col-1))+3:numRows
             C(row,col) = 1-shadeIncr*(row-(numRows-(col-1)));
        end

    end
    C = C*255;
    %imagesc(C)
    imagesc('XData',[1 6],'YData',[6 1],'CData',C) 

    %cmap = gray; 
    %cmap = cmap(20:64,:);
    %colormap(flip(cmap))
    colormap(cmap)
    ax = gca;
    
%% Lines 
    %Diagonal
    plot([0.5,6],[0.5,6],'--','LineWidth',2,'Color','k')
    
    %Line at 1 
    xlim = ax.XLim;
    plot([.5,6],[1,1],'--','LineWidth',2,'Color','k')
    
%% Axes
    ax.TickDir = 'out';
    ax.XLim = [0.5 6];
    ax.YLim = [0.5 6];
    ax.XTick = [1];
    ax.YTick = [1];



    