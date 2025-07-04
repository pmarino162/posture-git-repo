clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\kinematic comparison';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numIterations = 100; %Default - 10000
   cutoffNumTraj = 10; %Num trials that must be present in a condition to keep it for analysis 
   
%% Datasets to include in analysis 
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};         
    isoDatasetList = {'E20200116','E20200117','E20200120'};    

%% Main loop
    resultStruct = struct('animal',[],'dataset',[],'result',[]);
    structInd = 1;
    task = 'bci';
    for datasetList = {'N20171215','R20201020'}
        tic
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Subtract workspace centers
        [Data] = subtractReachWorkspaceCenters(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParamsKinematicComparison(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        
        % Only use N00 and A90 for Nigel to match Rocky's postures
        if strcmpi(dataset,'N20171215')
            trajStruct = trajStruct([trajStruct.posture] ~= 2);
        end
        
        % For earl reaching, keep only postures 1-4
        if ismember(dataset,{'E20210706','E20210707','E20210708'})
            trajStruct = trajStruct(ismember([trajStruct.posture],[1,2,3,4]));
        end

        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= cutoffNumTraj);    
        %get trajStructDims
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        %Get minimum number of condition trials and timestamps
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 
        % Get field name for comparisons below
        secondField = trajFields{2};
        capitalizedSecondField = [upper(secondField(1)) secondField(2:end)];
        fieldName = ['avg' capitalizedSecondField]; % The second field is the one we use, always making zSmoothFR the first field

        %% Preallocate sessionResultStruct - Compare posture A target A to posture B target B
        sessionResultStruct = struct('postureA',[],'targetA',[],'postureB',[],'targetB',[],'difference',[],'differenceAveraged',[]);
        sessionResultStructInd = 1;
        for postureA = postureList
            withinPostureTargetList = unique([trajStruct([trajStruct.posture]==postureA).target]);       
            for targetA = withinPostureTargetList
                for postureB = postureList(postureList>=postureA)
                    if any([trajStruct.posture]==postureB & [trajStruct.target]==targetA)
                        sessionResultStruct(sessionResultStructInd).postureA = postureA;
                        sessionResultStruct(sessionResultStructInd).targetA = targetA;
                        sessionResultStruct(sessionResultStructInd).postureB = postureB;
                        sessionResultStruct(sessionResultStructInd).targetB = targetA;
                        sessionResultStruct(sessionResultStructInd).difference = NaN(1,numIterations);
                        sessionResultStruct(sessionResultStructInd).differenceAveraged = NaN;
                        sessionResultStructInd = sessionResultStructInd + 1;
                    end
                end
            end
        end
        
        %% Do computation
        numPts = minNumTimestamps;
        for i = 1:numIterations
            %Split trajStruct into two groups
            [trajStruct1,trajStruct2] = splitDataforCVKinematicComparison(trajStruct,minNumCondTrials,binWidth); 
            %Make every comparison
            for postureA = postureList
                withinPostureTargetList = unique([trajStruct([trajStruct.posture]==postureA).target]);       
                for targetA = withinPostureTargetList
                    for postureB = postureList(postureList>=postureA)
                        if any([trajStruct.posture]==postureB & [trajStruct.target]==targetA)
                            
                            % Get trajectories for comparison
                            cond1 = find([trajStruct1.target]==targetA & [trajStruct1.posture]==postureA);
                            traj1 = trajStruct1(cond1).(fieldName).traj;  % Predicted trajectory
                            if postureA == postureB % Use trajStruct2 for within-condition comparisons
                                cond2 = find([trajStruct2.target]==targetA & [trajStruct2.posture]==postureB);
                                traj2 = trajStruct2(cond2).(fieldName).traj;  % Comparison trajectory
                            else % Use trajStruct1 for across-condition comparisons
                                cond2 = find([trajStruct1.target]==targetA & [trajStruct1.posture]==postureB);
                                traj2 = trajStruct1(cond2).(fieldName).traj;  % Comparison trajectory
                            end

                            % For the isometric force task, only use the
                            % first 3 columns (x,y,z) of force (not 4th
                            % col, which is total force)
                            if ismember(dataset,{'E20200116','E20200117','E20200120'})
                                traj1 = traj1(:,1:3);
                                traj2 = traj2(:,1:3);
                            end

                            % Get tempNumPts
                            if min([size(traj1,1),size(traj2,1)]) < numPts
                                tempNumPts = min([size(traj1,1),size(traj2,1)]);
                            else
                                tempNumPts = numPts;
                            end
                            %Compute comparison difference
                            difference = getMeanDist(traj1,traj2,tempNumPts);
                            %Store results in sessionResultStruct
                            sessionResultStructInd = find([sessionResultStruct.postureA]==postureA & ...
                                [sessionResultStruct.targetA]==targetA & [sessionResultStruct.postureB]==postureB);
                            sessionResultStruct(sessionResultStructInd).difference(i) = difference;                 
                        end
                    end
                end
            end
        end
        
        %Average across bootstraps
        for i = 1:length(sessionResultStruct)
            sessionResultStruct(i).differenceAveraged = mean([sessionResultStruct(i).difference]);
        end
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
        
        % Timing
        dataset
        numIterations
        toc
    end

%% Collect results for each monkey
    %Get monkeyList
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    %Create monkey result struct 
    monkeyResultStruct = struct('monkey',[],'acrossPostureDifference',[]);
    monkeyInd = 1;
    for monkey = monkeyList
      tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
      acrossPostureDifference = struct('postureDiff',[],'difference',[]);
      
      for session = 1:numel(tempResultStruct)
          dataset = tempResultStruct(session).dataset;
          result = tempResultStruct(session).result;
          postureA = [result.postureA];
          postureB = [result.postureB];
          targetA = [result.targetA];
          
          if ismember(dataset,{'N20171215','N20180221'}) %Nigel BCI: all postures are '1 away'
               postureRange = 1;
          else
              postureRange = max(postureA) - min(postureA);
          end
          
          for postureDiff = 0:1:postureRange
              
              % Collect trajectory differences
              if ismember(dataset,{'N20171215','N20180221'}) %Nigel BCI: all postures are '1 away'
                  difference = [result((abs(postureB-postureA)>0)==postureDiff).differenceAveraged];
              else
                  difference = [result(abs(postureB-postureA)==postureDiff).differenceAveraged];
              end
              
              % Store in struct
              if ismember(postureDiff, [acrossPostureDifference.postureDiff])
                  acrossPostureDifferenceIdx = find([acrossPostureDifference.postureDiff]==postureDiff);
                  acrossPostureDifference(acrossPostureDifferenceIdx).difference = [acrossPostureDifference(acrossPostureDifferenceIdx).difference,difference];
              else
                  if (size(acrossPostureDifference,2)==1) && (isempty(acrossPostureDifference(1).difference))
                      acrossPostureDifferenceIdx = 1;
                  else
                    acrossPostureDifferenceIdx = size(acrossPostureDifference,2) + 1;
                  end
                  acrossPostureDifference(acrossPostureDifferenceIdx).postureDiff = postureDiff;
                  acrossPostureDifference(acrossPostureDifferenceIdx).difference = difference;
              end
          end
      end
        monkeyResultStruct(monkeyInd).monkey = monkey;
        monkeyResultStruct(monkeyInd).acrossPostureDifference = acrossPostureDifference;
        monkeyInd = monkeyInd + 1;
    end
      
%% Run statistical tests
    for monkey = monkeyList
        acrossPostureDifference = monkeyResultStruct(strcmp([monkeyResultStruct.monkey],monkey{1,1})).acrossPostureDifference;
        postureDiffList = [acrossPostureDifference.postureDiff];
        postureDiffInd = 1;
        for postureDiff = postureDiffList(1:length(postureDiffList)-1)
            x = acrossPostureDifference([acrossPostureDifference.postureDiff]==postureDiff).difference;
            y = acrossPostureDifference([acrossPostureDifference.postureDiff]==(postureDiff+1)).difference;
            [h,p] = ttest2(x,y,'Tail','left');
            monkeyResultStruct(strcmp([monkeyResultStruct.monkey],monkey{1,1})).acrossPostureDifference(postureDiffInd).h_nextPosture = h;
            monkeyResultStruct(strcmp([monkeyResultStruct.monkey],monkey{1,1})).acrossPostureDifference(postureDiffInd).p_nextPosture = p;
            postureDiffInd = postureDiffInd + 1;
        end
        
    end
      
    if saveFig
        save(fullfile(saveDir,[task,'.mat']), "monkeyResultStruct");
    end
        

%% Plot histograms
    monkeyInd = 1;
    for monkey = monkeyList
        f=figure; hold on;
        figWidth = 175;
        figHeight = 150;
        if strcmpi(task, 'reach')
            figHeight = (4/5) * figHeight;
        end
        
        f.Position = [200 200 figWidth figHeight];
        fs = 5;
        acrossPostureDifference = monkeyResultStruct(strcmp([monkeyResultStruct.monkey],monkey{1,1})).acrossPostureDifference;
        postureDiffList = [acrossPostureDifference.postureDiff];
       
        for postureDiff = postureDiffList
           ax(postureDiff+1) = subplot(length(postureDiffList),1,postureDiff+1);
           difference = acrossPostureDifference([acrossPostureDifference.postureDiff]==postureDiff).difference;
           mean_difference = mean(difference);
           histogram(difference, 100, 'Normalization', 'probability');
           xline(mean_difference, 'r', 'LineWidth', 2);
           
           %title(['Posture difference: ',num2str(postureDiff)])
           %ylabel('Probability');
           if postureDiff ~= postureDiffList(end)
               xticklabels({});
           end
           set(gca,'fontname','arial')
           set(gca,'fontsize',fs)
        end
        if strcmpi(task, 'bci')
            %xlabel('Mean Euclidean distance between trajectories (mm)');
        elseif strcmpi(task, 'iso')
            %xlabel('Mean Euclidean distance between trajectories (N)');
        elseif strcmpi(task, 'reach')
            %xlabel('Mean Euclidean distance between trajectories (mm)');
        end
            
        linkaxes(ax, 'x');
        
        if saveFig
            saveas(gcf,fullfile(saveDir,['monkey_',monkey{1,1},'_',task,'.svg']));
            saveas(gcf,fullfile(saveDir,['monkey_',monkey{1,1},'_',task,'.png']));
        end
        
        clearvars ax;
    end


%% Local functions  
    %Get mean dist
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end
    
  
    