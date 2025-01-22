clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\Reviewer responses\Analyses\Fig 6\simulation_for_panel_b';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Generate trajectories
    numPts = 10;
    [traj1, traj2] = generateTrajectories();
    
%% Measure the difference before and after shift
    errorBeforeShift = getMeanDist(traj1,traj2,numPts);
    alpha = getOptimalAlpha(traj1,traj2,numPts);
    shiftTraj2 = traj2 + alpha;
    errorAfterShift = getMeanDist(traj1,shiftTraj2,numPts);
    
    
%% Plot trajectories before shift
    % Plot the trajectories excluding the last point
    figure;
    hold on;
    plot(traj1(1:end-1,1), traj1(1:end-1,2), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    plot(traj2(1:end-1,1), traj2(1:end-1,2), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
    
    % Plot the last segment without marker to avoid overlap with arrow
    plot(traj1(end-1:end,1), traj1(end-1:end,2), 'r-', 'LineWidth', 1.5);
    plot(traj2(end-1:end,1), traj2(end-1:end,2), 'b-', 'LineWidth', 1.5);

    addArrow(traj1, 'r');
    addArrow(traj2, 'b');
    
    axis equal;
    xlabel('Dim 1');
    ylabel('Dim 2');
    
    title(['\fontsize{16}Difference before translation: ', num2str(errorBeforeShift, '%.2f')]);
    legend('Trajectory 1', 'Trajectory 2 (before translation)');
    grid on;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,'before_translation_v2.jpg'));
    end

%% Plot trajectories after shift  
    % Plot the trajectories excluding the last point
    figure;
    hold on;
    plot(traj1(1:end-1,1), traj1(1:end-1,2), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    plot(shiftTraj2(1:end-1,1), shiftTraj2(1:end-1,2), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'b');

    % Plot the last segment without marker to avoid overlap with arrow
    plot(traj1(end-1:end,1), traj1(end-1:end,2), 'r-', 'LineWidth', 1.5);
    plot(shiftTraj2(end-1:end,1), shiftTraj2(end-1:end,2), 'b-', 'LineWidth', 1.5);

    addArrow(traj1, 'r');
    addArrow(shiftTraj2, 'b');
    axis equal;
    
    xlabel('Dim 1');
    ylabel('Dim 2');
    
    title(['\fontsize{16}Difference after translation: ', num2str(errorAfterShift, '%.2f')]);
    legend('Trajectory 1', 'Trajectory 2 (after translation)');
    grid on;


    if saveFig
        saveas(gcf,fullfile(saveDir,'after_translation_v2.jpg'));
    end
    
%% Local functions   
    % Create trajectories
   function [traj1, traj2] = generateTrajectories()
        % generateTrajectories generates two trajectories, traj1 and traj2.
        % Each trajectory is a 10x2 matrix where each row corresponds to a point on the trajectory.
        % Both trajectories start at (0,0) and move to the right, ending at x=10.
        % traj1 curves upwards, traj2 curves downwards.

        % Generate 10 x-coordinates from 0 to 10
        x = linspace(0, 10, 10);
        % Set amplitude for the curves
        amplitude = 4;

        % Calculate y-values for the upward-curving trajectory
        y1 = (x / 10).^3 * amplitude;
        % Calculate y-values for the downward-curving trajectory
        y2 = -(x / 10).^3 * amplitude;

        % Combine x and y-values into 10x2 matrices for each trajectory
        traj1 = [x', y1'];
        traj2 = [x', y2'];
    end

    % Add arrows to plot
    function addArrow(traj, colorStr)
        % Compute direction vectors at the last points
        dx = traj(end,1) - traj(end-1,1);
        dy = traj(end,2) - traj(end-1,2);
        dir = [dx, dy] / norm([dx, dy]);  % Normalize direction vector

        % Define arrow length
        arrowLength = norm(traj(end,:)-traj(end-1,:));

        % Compute arrow vectors
        dx_arrow = dir(1) * arrowLength;
        dy_arrow = dir(2) * arrowLength;

        % Plot arrows at the last point of each trajectory
        quiver(traj(end-1,1), traj(end-1,2), dx_arrow, dy_arrow, 'MaxHeadSize', 1.5, 'Color', colorStr, 'LineWidth', 1.5, 'AutoScale', 'off');
    end
    
    % Get optimal alpha
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end
    
    % Get mean dist
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end