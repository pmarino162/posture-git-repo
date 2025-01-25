function [moveOnsetTime,peakSpeedTime,rmTrialFlag] = getReachKin(trialData,trialInclStates,binWidth,kernelStdDev)    

%Gets time of movement onset and peakSpeed for center-out reaching style
%trials.  If the trial is ill-behaved, set rmTrialFlag to true to flag for
%subsequent removal.

    %Get marker velocity; peakSpeed; moveOnsetSpeed
    [markerVel,velTime] = getStatesTraj20220419(trialData,trialInclStates,'markerVel',binWidth,kernelStdDev,'timeRelToTrialStart',true);
    speed = vecnorm(markerVel');
    [peakSpeed,peakSpeedInd] = max(speed);
    peakSpeedTime = velTime(peakSpeedInd);
    moveOnsetSpeed = 0.2*peakSpeed;

    %Walk backwards from peak speed to get movement onset index and time
    rmTrialFlag = false;
    i = peakSpeedInd;
    while speed(i) > moveOnsetSpeed
        i = i-1;
        %remove trial if speed never dips below threshold during
        %reach epoch
        if i==0
            rmTrialFlag = true;
            i=1;
            break
        end
    end
    moveOnsetInd = i;
    moveOnsetTime = velTime(moveOnsetInd);  

    %Visualize trial velocity profile for debugging
%     [markerPos,posTime] = getStatesTraj20220419(trialData,trialInclStates,'markerPos',binWidth,kernelStdDev,'timeRelToTrialStart',true);
%                 %Test
%                 figure
%                     plot(velTime,markerVel(:,1))
%                     hold on
%                     plot(velTime,markerVel(:,2))
%                     ax = gca; xlim = ax.XLim; ylim = ax.YLim;
%                     title('velocity')
%                     legend('x','y')
%                     plot([moveOnsetTime moveOnsetTime],[ylim(1) ylim(2)],'--','Color','r')
% 
%                 figure
%                     plot(velTime,speed)
%                     hold on
%                     ax = gca; xlim = ax.XLim; ylim = ax.YLim;
%                     plot([xlim(1) xlim(2)],[moveOnsetSpeed moveOnsetSpeed],'--','Color','r')
%                     plot([moveOnsetTime moveOnsetTime],[ylim(1) ylim(2)],'--','Color','r')
%                     title('speed')
%                                  
%                     
%             figure; 
%             subplot(5,1,1)
%                 plot(velTime,speed); hold on;
%                 ax = gca; xlim = ax.XLim; ylim = ax.YLim;
%                 plot([xlim(1) xlim(2)],[moveOnsetSpeed moveOnsetSpeed],'--','Color','r')
%                 plot([moveOnsetTime moveOnsetTime],[ylim(1) ylim(2)],'--','Color','r')
%                 xticklabels();
%                 ylabel('speed (mm/ms)')
%             subplot(5,1,2)
%                 plot(velTime,markerVel(:,1)); hold on;
%                 ylabel('V_X (mm/ms)')
%             subplot(5,1,3)
%                 plot(velTime,markerVel(:,2)); hold on;
%                 ylabel('V_Y (mm/ms)')
%             subplot(5,1,4)
%                 plot(posTime,markerPos(:,1)); hold on;
%                 ylabel('x (mm)')
%             subplot(5,1,5)
%                 plot(posTime,markerPos(:,2)); hold on;
%                 ylabel('y (mm)')
%                 xlabel('time (ms)')
%         a = 1;
end