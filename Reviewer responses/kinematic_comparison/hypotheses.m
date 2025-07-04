clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\kinematic comparison';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numObs = 1000;
   var = 1;
   sigma = sqrt(var);
   x = linspace(-10, 10, 1000);
   postureDiffList = [0,1,2];
   
%% Create structs for each hypothesis
    hyp1Struct = struct('postureDiff',[],'u',[],'dist',[]);
    hyp2Struct = struct('postureDiff',[],'u',[],'dist',[]);
   
%% Fill structs
    % Hypothesis 1
    structInd = 1;
    for postureDiff = postureDiffList
        u = 1;
        dist = (1 / (sigma * sqrt(2*pi))) * exp( - (x - u).^2 ./ (2*var) );
        hyp1Struct(structInd).postureDiff = postureDiff;
        hyp1Struct(structInd).u = u;
        hyp1Struct(structInd).dist = dist;
        structInd = structInd + 1;
    end

    %Hypothesis 2
    structInd = 1;
    for postureDiff = postureDiffList
        u = 1+ 2*postureDiff;
        dist = (1 / (sigma * sqrt(2*pi))) * exp( - (x - u).^2 ./ (2*var) );
        hyp2Struct(structInd).postureDiff = postureDiff;
        hyp2Struct(structInd).u = u;
        hyp2Struct(structInd).dist = dist;
        structInd = structInd + 1;
    end
      
%% Plot histograms
% Hyp 1
    fig=figure; hold on;
        figWidth = 175;
        figHeight = 150;
        fig.Position = [200 200 figWidth figHeight];
    for postureDiff = postureDiffList
        ax(postureDiff+1) = subplot(length(postureDiffList),1,postureDiff+1);
        dist = hyp1Struct([hyp1Struct.postureDiff]==postureDiff).dist;
        u = hyp1Struct([hyp1Struct.postureDiff]==postureDiff).u;
        plot(x, dist, 'LineWidth', 2);
        hold on;
        fill([x, fliplr(x)], [dist, zeros(size(dist))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xline(u, 'r', 'LineWidth', 3);
        xticks([]);
        yticks([]);
        xlim([-3, 9]);
    end
    linkaxes(ax, 'x');

    if saveFig
        saveas(gcf,fullfile(saveDir,['hyp1.svg']));
        saveas(gcf,fullfile(saveDir,['hyp1.png']));
    end
   
% Hyp 2
        
    fig=figure; hold on;
        figWidth = 175;
        figHeight = 150;
        fig.Position = [200 200 figWidth figHeight];
    for postureDiff = postureDiffList
        ax(postureDiff+1) = subplot(length(postureDiffList),1,postureDiff+1);
        dist = hyp2Struct([hyp2Struct.postureDiff]==postureDiff).dist;
        u = hyp2Struct([hyp2Struct.postureDiff]==postureDiff).u;
        plot(x, dist, 'LineWidth', 2);
        hold on;
        fill([x, fliplr(x)], [dist, zeros(size(dist))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xline(u, 'r', 'LineWidth', 3);
        xticks([]);
        yticks([]);
        xlim([-3, 9]);
    end
    linkaxes(ax, 'x');
    
    if saveFig
        saveas(gcf,fullfile(saveDir,['hyp2.svg']));
        saveas(gcf,fullfile(saveDir,['hyp2.png']));
    end
