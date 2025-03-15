clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = "C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\kinematic comparison\trajectory_plots";
    set(0, 'DefaultFigureRenderer', 'painters');
   
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};   
    
%% Setup colormap (based on number of postures)
    [pcmap,tcmap,rainbow] = getColorMaps(5); 
    
%% Load data, getTrajStruct
    dataset = 'E20210706';
    monkey = dataset(1);
    task = 'reach';
    %Load data
    [Data,zScoreParams] = loadData(dataset);
    [Data] = removeShortBCIandIsoTrials(Data,dataset);
    %Subtract workspace centers
    [Data] = subtractReachWorkspaceCenters(Data,dataset);
    
    %Get trajStruct
    [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParamsKinematicComparison(dataset);
    trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
    %Remove any conditions for which there weren't enough trials
    [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);            
    %get trajStructDims
    [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
    %Get minimum number of condition trials and timestamps
    [minNumCondTrials] = getMinNumCondTrials(trajStruct);
    [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
    % Get field name for comparisons below
    secondField = trajFields{2};
    capitalizedSecondField = [upper(secondField(1)) secondField(2:end)];
    fieldName = ['avg' capitalizedSecondField]; % The second field is the one we use, always making zSmoothFR the first field
    allFieldName = ['all' capitalizedSecondField];

    
%% Plot Positions using trajStruct 
fs = 5;
    for posture = postureList

        fig = figure; hold on
        fig.Position = [200 200 100 100];
        
        for target = targetList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target);
            if ~isempty(tempTrajStruct)

                % Plot individual trial trajectories (faded)
                numTrials = size(tempTrajStruct.(allFieldName),2);
                for trial = 1:numTrials
                    traj = tempTrajStruct.(allFieldName)(trial).traj;
                    %plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
                    plot(traj(:,1), traj(:,2), ...
                     'Color', [tcmap(target,:) 0.2], ... 
                     'LineWidth', 0.25);
                end
                
                % Plot average trajectory (bold)                
                traj = tempTrajStruct.(fieldName).traj;                
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:), 'LineWidth', 0.5);
                
            end
        end
        

        %axis('square');
        axis('equal');
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        
        
        title(['Posture ',num2str(posture)])
        xlabel('mm');
        ylabel('mm');
        if saveFig
            saveas(gcf,fullfile(saveDir,['monkey_',monkey,'_',task,'posture_',num2str(posture),'.svg']));
            saveas(gcf,fullfile(saveDir,['monkey_',monkey,'_',task,'posture_',num2str(posture),'.png']));
        end
    end
    
    
    
%     for posture = postureList
%         fig = figure; hold on
%         fig.Position = [200 200 50 50];
%         for target = targetList
%             tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target);
%             if ~isempty(tempTrajStruct)
% 
%                 % Plot individual trial trajectories (faded)
%                 numTrials = size(tempTrajStruct.(allFieldName),2);
%                 for trial = 1:numTrials
%                     traj = tempTrajStruct.(allFieldName)(trial).traj;
%                     time = tempTrajStruct.(allFieldName)(trial).timestamps;
%                     for dim = 1:3
%                         subplot(3,1,dim)
%                         hold on
%                         plot(time, traj(:,dim),'Color', [tcmap(target-2,:) 0.2],'LineWidth', 0.25);
%                     end
%                 end
%                 
%                 % Plot average trajectory (bold)                
%                 traj = tempTrajStruct.(fieldName).traj; 
%                 time = tempTrajStruct.(fieldName).timestamps;
%                 for dim = 1:3
%                     subplot(3,1,dim)
%                     hold on
%                     plot(time, traj(:,dim),'Color',tcmap(target-2,:), 'LineWidth', 0.5);
%                     if dim == 1
%                         title(['Posture ',num2str(posture)])
%                         %ylabel('Force X (N)');
%                         
%                         xticks([]);
%                     elseif dim == 2
%                         %ylabel('Force Y (N)');
%                         xticks([]);
%                     elseif dim == 3
%                         %ylabel('Force Z (N)');
%                          %xlabel('time (ms)');
%                     end
%                     
%                     ylims = ylim;
% 
%                     % Create tick locations from the actual min and max
%                     yTicks = [ylims(1), ylims(2)];
% 
%                     % Set those ticks, but round their labels to the nearest whole number
%                     set(gca, 'YTick', yTicks, ...
%                              'YTickLabel', {num2str(round(ylims(1))), num2str(round(ylims(2)))});
% 
%                     % (Optional) turn off minor ticks
%                     set(gca, 'YMinorTick','off');
%                     
%                     set(gca,'fontname','arial')
%                     set(gca,'fontsize',fs)
%                 end
%                 
% 
%                 
%             end
%         end
%         
%         
%         if saveFig
%             saveas(gcf,fullfile(saveDir,['monkey_',monkey,'_',task,'posture_',num2str(posture),'.svg']));
%             saveas(gcf,fullfile(saveDir,['monkey_',monkey,'_',task,'posture_',num2str(posture),'.png']));
%         end
%     end