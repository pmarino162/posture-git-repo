clear; clc; clf; close all

%% Load dataset
    dataset = 'E20210706';
    [Data,zScoreParams] = loadData(dataset);

%% Histogram of rxnTimes
    kinData = [Data.kinData];
    rxnTime = [kinData.rxnTime];
    reachTime = [kinData.reachTime];
    
%% Identify trials with weird rxn times & reach times 
    figure; histogram(rxnTime); xlabel('rxn time (ms)'); ylabel('count');
    figure; histogram(reachTime); xlabel('reach time (ms)'); ylabel('count');
    
%% Visualize position, velocity, and speed for that trial

    trialsOfInterest = find(rxnTime > 100 & rxnTime < 150);
    
    %trial = 400;
    trial = 43;
    
    %Get marker velocity; peakSpeed; moveOnsetSpeed
    binWidth = 25; kernelStdDev = 25; trajFields = {'markerPos','markerVel'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        case{'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};  
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Target Acquire','first',0},{'state','Target Hold','first',0}};
        case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
            trialInclStates(1).trialName = {'Nigel Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Target Hold','first',0}};
        case{'R20200221','R20200222'}
            trialInclStates(1).trialName = {'Rocky Dissociation'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Target Hold','first',0}};
    end
    [markerVel,velTime] = getStatesTraj20220419(Data(trial),trialInclStates,'markerVel',binWidth,kernelStdDev,'timeRelToTrialStart',true);
    speed = vecnorm(markerVel');
    [peakSpeed,peakSpeedInd] = max(speed);
    peakSpeedTime = velTime(peakSpeedInd);
    moveOnsetSpeed = 0.2*peakSpeed;
      
    moveOnsetTime = Data(trial).kinData.moveOnsetTime;
    rxnTime = Data(trial).kinData.rxnTime
    %Visualize trial velocity profile for debugging
        %Test
        figure
            plot(velTime,markerVel(:,1))
            hold on
            plot(velTime,markerVel(:,2))
            ax = gca; xlim = ax.XLim; ylim = ax.YLim;
            title('velocity')
            legend('x','y')
            plot([moveOnsetTime moveOnsetTime],[ylim(1) ylim(2)],'--','Color','r')
        figure
            plot(velTime,speed)
            hold on
            ax = gca; xlim = ax.XLim; ylim = ax.YLim;
            plot([xlim(1) xlim(2)],[moveOnsetSpeed moveOnsetSpeed],'--','Color','r')
            plot([moveOnsetTime moveOnsetTime],[ylim(1) ylim(2)],'--','Color','r')
            title('speed')