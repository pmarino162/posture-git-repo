function [Data] = applyExclusionCriteria(Data)

%Divide dataset up by task. Apply exclusion criteria to each task
%separately, then put back together in order. 

trialNames = {Data.trialName};

isoTrials = [cellfun(@(x) strcmpi(x,'IsometricForce_1D'),trialNames)]==1;

bciTrials = [cellfun(@(x) strcmpi(x,'GridTask_CO_Across_BC_ForceBar') || strcmpi(x,'GridTask_BC_ForceBar') ...
    ||strcmpi(x,'Nigel Posture BC Center Out') || strcmpi(x,'Rocky Posture BC Center Out') || strcmpi(x,'BCI Center Out') || strcmpi(x,'CenterOut_ForceBar_BC'), ...
    trialNames)]==1;

reachTrials = [cellfun(@(x) strcmpi(x,'GridReaching') || strcmpi(x,'HC_CenterOut_ForceBar_20200314') ...
    ||strcmpi(x,'Rocky Dissociation') || strcmpi(x,'Nigel Dissociation') ,trialNames)]==1;

bciData = Data(bciTrials);
isoData = Data(isoTrials);
reachData = Data(reachTrials);

%BCI : remove trials with moveTime >= 95th %ile
if ~isempty(bciData)
    kinData = [bciData.kinData];
    moveTime = [kinData.moveTime];
    moveTimeCutoff = prctile(moveTime,95);
    bciData = bciData(moveTime < moveTimeCutoff);
    %Histogram for validation
    figure; hold on;
        histogram(moveTime)
        ax = gca; ylimits = ax.YLim;
        plot(moveTimeCutoff*ones(1,2),ylimits,'-r');
        xlabel('Move Time'); ylabel('Number of trials')
end

%Isometric force : remove trials with reachTime >= 95th %ile
if ~isempty(isoData)
    kinData = [isoData.kinData];
    reachTime = [kinData.reachTime];
    reachTimeCutoff = prctile(reachTime,95);
    isoData = isoData(reachTime < reachTimeCutoff);
    %Histogram for validation
    figure; hold on;
        histogram(reachTime)
        ax = gca; ylimits = ax.YLim;
        plot(reachTimeCutoff*ones(1,2),ylimits,'-r');
        xlabel('Reach Time'); ylabel('Number of trials')
end

%Multiple tasks reaching - 
%Reaching : remove trials with reachTime > 95th %ile
if ~isempty(reachData)
    kinData = [reachData.kinData];
    reachTime = [kinData.reachTime];
    reachTimeCutoff = prctile(reachTime,95);
    reachData = reachData(reachTime < reachTimeCutoff);
    %Histogram for validation
    figure; hold on;
        histogram(reachTime)
        ax = gca; ylimits = ax.YLim;
        plot(reachTimeCutoff*ones(1,2),ylimits,'-r');
        xlabel('Reach Time'); ylabel('Number of trials')
end

%Put back together and sort by trial number
Data = [bciData,isoData,reachData];
trialNum = [Data.trialNum];
[~,sortInd] = sort(trialNum);
Data = Data(sortInd);

end