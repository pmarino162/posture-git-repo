function [moveOnsetTime,rmTrialFlag] = getIsoKin(trialData,binWidth,kernelStdDev)    

%Gets time of movement onset for isometric force style
%trials.  If the trial is ill-behaved, set rmTrialFlag to true to flag for
%subsequent removal.

    %Get trialName info for trialInclStates
    trialName = trialData.trialName;
    trialInclStates = struct('trialName','','inclStates',[]);
    trialInclStates(1).trialName = {trialName};
    
    %Get resting force
    trialInclStates(1).inclStates = {{'state','Center Hold','last',25},{'state','Target','first',25}};
    [restForce,restForceTime] = getStatesTraj20220419(trialData,trialInclStates,'force',binWidth,kernelStdDev,'timeRelToTrialStart',true);

    %Get movement force
    trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
    [force,forceTime] = getStatesTraj20220419(trialData,trialInclStates,'force',binWidth,kernelStdDev,'timeRelToTrialStart',true);
    
    %Get "yForce" = abs(yForceMove-yForceRest)
    yForce = abs(force(:,2)-mean(restForce(:,2)));
    
    %Get y-force; max y-force; moveOnsetForce
    [peakForce,peakForceInd] = max(yForce);
    peakForceTime = forceTime(peakForceInd);
    moveOnsetForce = 0.2*peakForce;

    %Walk backwards from peak speed to get movement onset index and time
    rmTrialFlag = false;
    i = peakForceInd;
    while yForce(i) > moveOnsetForce
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
    moveOnsetTime = forceTime(moveOnsetInd);  

%     %Visualize trial force profile for debugging
%             %Test
%             figure; 
%             subplot(5,1,1)
%                 plot(forceTime,yForce); hold on;
%                 ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
%                 plot(xlimits,moveOnsetForce*ones(1,2),'-- r')
%                 plot(moveOnsetTime*ones(1,2),ylimits,'-r')
%                 xticklabels();
%                 ylabel('yForce (N)')
%             subplot(5,1,2)
%                 plot(forceTime,force(:,4)); hold on;
%                 ylabel('Total Force (N)')
%             subplot(5,1,3)
%                 plot(forceTime,force(:,1)); hold on;
%                 ylabel('F_X (N)')
%             subplot(5,1,4)
%                 plot(forceTime,force(:,2)); hold on;
%                 ylabel('F_Y (N)')
%             subplot(5,1,5)
%                 plot(forceTime,force(:,3));hold on; 
%                 ylabel('F_Z (N)')
%                 xlabel('time (ms)')
%                 
%                 a=1;
end