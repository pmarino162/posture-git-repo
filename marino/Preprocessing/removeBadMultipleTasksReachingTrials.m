function [Data] = removeBadMultipleTasksReachingTrials(Data)    

%Removes trials for which hand was not resting on forcebar before reach
%began

%% Identify ill-behaved trials
numTrials = size(Data,2);
binWidth = 1; %Necessary input to getStatesTraj20220419, but not used
kernelStdDev = NaN; %Necessary input to getStatesTraj20220419, but not used
rmTrials = [];

%Set threshold for removal (N)
hiThreshold = 5;
loThreshold = 1.5;

%Set up parts of trialInclStates which will be used throughout
trialInclStates = struct('trialName','','inclStates',[]);
trialInclStates(1).trialName = {'HC_CenterOut_ForceBar_20200314'};
trialInclStates(1).inclStates = {{'state','Reach','first',-200},{'state','Reach','first',0}};

for trial = 1:numTrials
    %Only conisder reaching trials
    trialData = Data(trial);
    trialName = trialData.trialName;
    
    if strcmpi(trialName,'HC_CenterOut_ForceBar_20200314')
        %Get force around go cue
        [force,forceTime] = getStatesTraj20220419(trialData,trialInclStates,'force',binWidth,kernelStdDev,'timeRelToTrialStart',true);

        %Get mean force magnitdue (away from desired force of 2N)
        meanForceMag = mean(force(:,4));
        
        %If force exceeds threshold, flag trial for removal
        if meanForceMag > hiThreshold || meanForceMag < loThreshold
           rmTrials = [rmTrials,trial]; 
        end

        %Visualize forces (optional)
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
    
end

%% Remove ill-behaved trials
    Data(rmTrials) = [];
    


   


end