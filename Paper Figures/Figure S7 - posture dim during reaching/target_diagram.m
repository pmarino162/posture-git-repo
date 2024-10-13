clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure S7 - posture dim during reaching';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Parameters
    targetMarkerSize = 80;
    
%% Monkey E
    workspaceCenter1 = [-40.0, 40.0];
    workspaceCenter2 = [40.0, -40.0];
    targetOffset = sqrt((80^2)/2);
    targetOffsetVec = [targetOffset, -targetOffset];

    targetPositions = [ workspaceCenter1-targetOffsetVec;
                        workspaceCenter1; ...
                        workspaceCenter1+targetOffsetVec; ...
                        workspaceCenter2-targetOffsetVec;
                        workspaceCenter2; ...
                        workspaceCenter2+targetOffsetVec];
 %% Plot
    figure; hold on;
    
    % Targets
    numTargets = size(targetPositions,1);
    for i=1:numTargets
        plot(targetPositions(i,1), targetPositions(i,2), '.', 'MarkerSize', targetMarkerSize, 'Color', 'k')
    end
    
    % Axis
    offset1 = 40;
    plot([targetPositions(1,1)-targetOffsetVec(1)/2, targetPositions(end,1)+targetOffsetVec(1)/2]+offset1,...
        [targetPositions(1,2)-targetOffsetVec(2)/2, targetPositions(end,2)+targetOffsetVec(2)/2]+offset1,'LineWidth',2, 'Color','k')
    
    % Ticks
    offset2 = 20;
    for i = 1:6
        tempPos = [targetPositions(i,1), targetPositions(i,2)];
        plot([tempPos(1), tempPos(1)+offset2]+offset1,[tempPos(2), tempPos(2)+offset2]+offset1,'LineWidth',2, 'Color','k')
    end
    plot([0, offset2]+offset1,[0, offset2]+offset1,'LineWidth',2, 'Color','k')
    axis equal
    if saveFig
        saveas(gcf,fullfile(saveDir,'Monkey_E_diagram.svg'));
    end                
                    
%% Monkeys N and R
%% Get target positions
    workspaceCenter1 = [-51.2, 34.7];
    workspaceCenter2 = [51.2, -34.7];
    targetOffset = sqrt((60^2)/2);
    targetOffsetVec = [targetOffset, -targetOffset];

    targetPositions = [ workspaceCenter1-targetOffsetVec;
                        workspaceCenter1; ...
                        workspaceCenter1+targetOffsetVec; ...
                        workspaceCenter2-targetOffsetVec;
                        workspaceCenter2; ...
                        workspaceCenter2+targetOffsetVec];

%% Plot
    figure; hold on;
    % Targets
    numTargets = size(targetPositions,1);
    for i=1:numTargets
        plot(targetPositions(i,1), targetPositions(i,2), '.', 'MarkerSize', targetMarkerSize, 'Color', 'k')
    end
    
    % Axis
    offset1 = 40;
    plot([targetPositions(1,1)-targetOffsetVec(1)/2, targetPositions(1,1)+170+targetOffsetVec(1)/2]+offset1,...
        [targetPositions(1,2)-targetOffsetVec(2)/2, targetPositions(1,2)-170+targetOffsetVec(2)/2]+offset1,'LineWidth',2, 'Color','k')
    
    % Ticks
    offset2 = 20;
    for i = 1:5
        tempPos = [targetPositions(1,1), targetPositions(1,2)] + (i-1)*targetOffsetVec;
        plot([tempPos(1), tempPos(1)+offset2]+offset1,[tempPos(2), tempPos(2)+offset2]+offset1,'LineWidth',2, 'Color','k')
    end
    axis equal     
    if saveFig
        saveas(gcf,fullfile(saveDir,'Monkeys_NR_diagram.svg'));
    end    