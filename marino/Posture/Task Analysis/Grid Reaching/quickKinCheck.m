function [] = quickKinCheck(trialData)
    
    %Set up trial include states
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {trialData.trialName};
    trialInclStates(1).addTimeToBeginning = {0};
    trialInclStates(1).addTimeToEnd = {0};
    trialInclStates(1).inclStates = {'Target Acquire'};
    trialInclStates(1).inclOccurrence = {'first'};
    
    %Get pos vel speed 
    [markerPos,posTime,~] = getStatesTraj(trialData,trialInclStates,'marker','timeRelToTrialStart',true);
    [markerVel,velTime,~] = getStatesTraj(trialData,trialInclStates,'markerVel','timeRelToTrialStart',true);
    speed = vecnorm(markerVel');
    [peakSpeed,peakSpeedInd] = max(speed);
    peakSpeedTime = velTime(peakSpeedInd);
    moveOnsetSpeed = 0.2*peakSpeed;
    %Walk backwards from peak speed to get movement onset index
    i = peakSpeedInd;
    while speed(i) > moveOnsetSpeed
        i = i-1;
    end
    moveOnsetInd = i;
    moveOnsetTime = velTime(moveOnsetInd);  

    %Plot 
    
    figure
    plot(velTime,markerPos(:,1))
    hold on
    plot(velTime,markerPos(:,2))
    title('position')
    
    
    figure
    plot(velTime,markerVel(:,1))
    hold on
    plot(velTime,markerVel(:,2))
    title('velocity')

    figure
    plot(velTime,speed)
    title('speed')




end