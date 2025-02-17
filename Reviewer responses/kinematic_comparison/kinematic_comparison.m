clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 6';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Set parameters
   numIterations = 20;
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
    for datasetList = {'E20200316','R20200221'}%{'R20200221','N20190222'}%{'N20190222','N20190226','R20200221','R20200222'}%reachDatasetList%{'E20200316'}%bciDatasetList% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %% Set up trajStruct
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        [Data] = removeShortBCIandIsoTrials(Data,dataset);
        %Get trajStruct
        [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParamsKinematicComparison(dataset);
        trajStruct = getTrajStruct(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'zScoreParams',zScoreParams,'getTrialAverages',false);  
        %Remove any conditions for which there weren't enough trials
        [numCondTrials] = getNumCondTrials(trajStruct,'showHist',false);         
        trajStruct = trajStruct(numCondTrials >= cutoffNumTraj);    
        %get trajStructDims
        [postureList,numPostures,targetList,numTargets,numChannels,numConditions] = getTrajStructDimensions(trajStruct);
        %Get minimum number of condition trials and timestamps
        [minNumCondTrials] = getMinNumCondTrials(trajStruct);
        [minNumTimestamps] = getMinNumTimestamps(trajStruct); 

                
        %% Preallocate sessionResultStruct - Compare posture A target A to posture B target B
        sessionResultStruct = struct('postureA',[],'targetA',[],'postureB',[],'targetB',[],'difference',[]);
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
                            cond1 = find([trajStruct1.target]==targetA & [trajStruct1.posture]==postureA);
                            cond2 = find([trajStruct1.target]==targetA & [trajStruct1.posture]==postureB);
                            %Compute comparison difference 
                            secondField = trajFields{2};
                            capitalizedSecondField = [upper(secondField(1)) secondField(2:end)];
                            fieldName = ['avg' capitalizedSecondField]; % The second field is the one we use, always making zSmoothFR the first field
                            
                            traj1 = trajStruct1(cond1).(fieldName).traj;  % Predicted trajectory
                            if postureA == postureB % Use trajStruct2 for within-condition comparisons
                                traj2 = trajStruct2(cond2).(fieldName).traj;  % Comparison trajectory
                            else
                                traj2 = trajStruct1(cond2).(fieldName).traj;  % Comparison trajectory
                            end


                            % Get tempNumPts
                            if min([size(traj1,1),size(traj2,1)]) < numPts
                                tempNumPts = min([size(traj1,1),size(traj2,1)]);
                            else
                                tempNumPts = numPts;
                            end
                            %Compute absolute and normalized errors
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

        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).result = sessionResultStruct;
        structInd = structInd + 1;
    end

%% Collect results for each monkey
    %Get monkeyList
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    %Create monkey result struct (reference posture 1)
    monkeyResultStruct = struct('monkey',[],'acrossPostureError',[]);
    monkeyInd = 1;
    for monkey = monkeyList
      tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
      acrossPostureDifference = struct('postureDiff',[],'difference',[]);
      acrossPostureDifferenceInd = 1;
      
      
      
      postureRange = 
      for postureDiff = 0:1:4
          difference = [];
          
          for session = 1:numel(tempResultStruct)
              result = tempResultStruct(session).result;
              postureA = [result.postureA];
              postureB = [result.postureB];
              targetA = [result.targetA];
              
              if postureDiff == 0
                 % Get first instance of each unique combo
                 [uniquePairs, keepIdx] = unique([postureA',targetA'], 'rows', 'first');
                 temp = [result(keepIdx).baselineError];
              else
                 % Include only comparisons to reference posture
                 temp = [result(abs(postureB-postureA)==postureDiff).difference];
              end
              
              difference = [difference, temp];
          end
          
          
          acrossPostureDifference(acrossPostureDifferenceInd).postureDiff = postureDiff;
          acrossPostureDifference(acrossPostureDifferenceInd).difference = difference;
          acrossPostureDifferenceInd = acrossPostureDifferenceInd + 1;
      end
      
      

      monkeyResultStruct(monkeyInd).monkey = monkey;
      monkeyResultStruct(monkeyInd).acrossPostureError = acrossPostureDifference;
      monkeyInd = monkeyInd + 1;
    end
    
%% Plot histograms
    monkeyInd = 1;
    for monkey = monkeyList
        f=figure; hold on;
        %monkeyResult = monkeyResultStruct(strcmp(monkeyResultStruct.monkey,monkey{1,1}).result;
        acrossPostureDifference = monkeyResultStruct.acrossPostureError;
        for postureDiff = 0:1:4
           subplot(5,1,postureDiff+1)
           error = acrossPostureDifference([acrossPostureDifference.postureDiff]==postureDiff).error;
           histogram(error, 'Normalization', 'probability');
        end
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
    
  
    