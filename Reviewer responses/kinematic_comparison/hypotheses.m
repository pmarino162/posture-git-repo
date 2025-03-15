clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\kinematic comparison';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numObs = 1000;
   sigma = 1;
   postureDiffList = [0,1,2];
   
%% Create structs for each hypothesis
    hyp1Struct = struct('postureDiff',[],'dist',[]);
    hyp2Struct = struct('postureDiff',[],'dist',[]);
   
%% Fill structs
    % Hypothesis 1
    structInd = 1;
    for postureDiff = postureDiffList
        mu = 1;
        dist = normrnd(mu, sigma, 1, numObs);
        hyp1Struct(structInd).postureDiff = postureDiff;
        hyp1Struct(structInd).dist = dist;
        structInd = structInd + 1;
    end

    %Hypothesis 2
    structInd = 1;
    for postureDiff = postureDiffList
        mu = 1+ 2*postureDiff;
        dist = normrnd(mu, sigma, 1, numObs);
        hyp2Struct(structInd).postureDiff = postureDiff;
        hyp2Struct(structInd).dist = dist;
        structInd = structInd + 1;
    end
      
%% Plot histograms
    fig=figure; hold on;
        figWidth = 175;
        figHeight = 150;
        fig.Position = [200 200 figWidth figHeight];
    for postureDiff = postureDiffList
       ax(postureDiff+1) = subplot(length(postureDiffList),1,postureDiff+1);
       difference = hyp1Struct([hyp1Struct.postureDiff]==postureDiff).dist;
       mean_difference = mean(difference);
        [f, x] = ksdensity(difference);
        plot(x, f, 'LineWidth', 2);
        hold on;
        fill([x, fliplr(x)], [f, zeros(size(f))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xline(mean_difference, 'r', 'LineWidth', 3);
        xticks([]);
        yticks([]);
        xlim([-4, 9]);
    end
%linkaxes(ax, 'x');

        if saveFig
            saveas(gcf,fullfile(saveDir,['hyp1.svg']));
            saveas(gcf,fullfile(saveDir,['hyp1.png']));
        end
        
    fig=figure; hold on;
            figWidth = 175;
        figHeight = 150;
        fig.Position = [200 200 figWidth figHeight];
    for postureDiff = postureDiffList
       ax(postureDiff+1) = subplot(length(postureDiffList),1,postureDiff+1);
       difference = hyp2Struct([hyp2Struct.postureDiff]==postureDiff).dist;
       mean_difference = mean(difference);
        [f, x] = ksdensity(difference);
        plot(x, f, 'LineWidth', 2);
        hold on;
        fill([x, fliplr(x)], [f, zeros(size(f))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xline(mean_difference, 'r', 'LineWidth', 3);
        xticks([]);
        yticks([]);
        xlim([-4, 9]);
    end
    %linkaxes(ax, 'x');
    
        if saveFig
            saveas(gcf,fullfile(saveDir,['hyp2.svg']));
            saveas(gcf,fullfile(saveDir,['hyp2.png']));
        end
