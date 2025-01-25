clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20230203\Supplemental figs\Hand Speed';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    tcmap = customRainbow;
    mcmap = tcmap([1,4,6],:);

%% Create allTrajStruct
    allTrajStruct = struct('dataset',[],'task',[],'trajStruct',[]);
    allTrajStructInd = 1;
    
%% Collect Monkey allTrajStructs
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    
    monkeyEDatasetList = {'E20200316','E20200317','E20200318','E20200319','E20210706','E20210707','E20210708'};
        monkey = 'E';
    
    monkeyNDatasetList = {'N20190222','N20190226','N20190227','N20190228','N20190307','N20171215','N20180221'};
        monkey = 'N';
    
    monkeyRDatasetList = {'R20200221','R20200222','R20201020','R20201021'};
        monkey = 'R';
    
    for datasetList = monkeyNDatasetList
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR','markerPos','markerVel'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316PSP_Sync','E20200316','E20200317','E20200318','E20200319'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
                task = 'bci';
                
                  case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';
                
                
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',300}};
                task = 'reach';      
                
                 case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
                
                
        end
        
        
        %Remove any long or short trials - Reach
        if strcmpi(task,'reach')
            kinData = [Data.kinData];
            rxnTime = [kinData.rxnTime];
            reachTime = [kinData.reachTime];
            figure
                histogram(rxnTime)
                xlabel('Reaction Time')
                ylabel('Number of conditions')
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of conditions')
           rxnCutoff = prctile(rxnTime,95);
           reachCutoff = prctile(reachTime,95);
           rmTrials =  unique([find(rxnTime > rxnCutoff), find(reachTime > reachCutoff)]);
           Data(rmTrials) = [];
       end

       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
        %Remove long/short trials, BCI
        if strcmpi(task,'bci')
            reachTime = [];
            for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
            end     
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of trials')
           reachCutoff = prctile(reachTime,95);
           for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               reachTime = [];
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
               rmTrials =  find(reachTime > reachCutoff);
               trajStruct(i).allSmoothFR(rmTrials) = [];
               trajStruct(i).allMarkerVel(rmTrials) = [];
            end  
        end

        allTrajStruct(allTrajStructInd).dataset = dataset;
        allTrajStruct(allTrajStructInd).task = task;
        allTrajStruct(allTrajStructInd).trajStruct = trajStruct;
        allTrajStructInd = allTrajStructInd + 1;
    end
    
%% Do Plotting 
for task = {'reach','bci'}
   %Combine all task data
        trajStruct = [];
        for i = 1:size(allTrajStruct,2)
            if strcmpi(allTrajStruct(i).task,task{1,1})
                trajStruct = [trajStruct,allTrajStruct(i).trajStruct];
            end
        end
   %Compute plot info
        allTraj = NaN(5000,100);
        allTrajInd = 1;
        trajLengths = [];
        startTimes = nan(1,5000);
        for i=1:size(trajStruct,2)
            numTraj = size(trajStruct(i).allMarkerVel,2);
            reachTime = [];
            for j = 1:numTraj
                traj = trajStruct(i).allMarkerVel(j).traj;
                timestamps = trajStruct(i).allMarkerVel(j).timestamps;
                if ~isempty(traj)
                    traj = vecnorm(traj')';
                    trajLength = size(traj,1);
                    trajLengths = [trajLengths,trajLength];
                    allTraj(allTrajInd,1:trajLength) = traj';
                    startTimes(1,allTrajInd) = timestamps(1);
                    allTrajInd = allTrajInd + 1;
                end
            end
        end
                
        medianTrajLength = median(trajLength);
        avgTraj = nanmean(allTraj);
        stdDev = nanstd(allTraj); 
        rtNumTrials = zeros(1,size(allTraj,2));
        for i = 1:size(allTraj,2)
            rtNumTrials(i) = sqrt(size(allTraj,1)-sum(isnan(allTraj(:,i))));
        end       
        stdError = stdDev./rtNumTrials;
        
        minTime = min(startTimes,[],'omitnan');
        time = minTime:binWidth:binWidth*(medianTrajLength-1);
   %Plot
        fs = 14;
        figure
        shadedErrorBar(time(1:medianTrajLength),avgTraj(1:medianTrajLength),stdError(1:medianTrajLength));
        xlabel('time after go cue (ms)')
        ylabel('hand speed (m/s)')
        set(gca, 'TickDir', 'out')
        set(gca,'fontname','arial'); set(gca,'fontsize',fs)
        
   %Save y lims if reach
        if strcmpi(task{1,1},'reach')
            ax = gca;
            reachYLims = ax.YLim;
        end
   %Use y lims if BCI
        if strcmpi(task{1,1},'bci')
            ylim(reachYLims);
        end
   %Save figure
        if saveFig
            saveas(gcf,fullfile(saveDir,['monkey',monkey,'_',task{1,1},'_handSpeed.svg']));
        end
end
        

