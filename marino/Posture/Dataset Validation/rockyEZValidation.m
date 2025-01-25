%% Load dataset
        EZData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390.raw_EZ.mat');
        EZ = EZData.EZ;
        
%% Unpack trial info
        clearvars -except EZ
        trial = 137;
        
        trialHeaderName = ['trial_header_',num2str(trial)];
        trial_header = EZ.(trialHeaderName);
        spikeDataName = ['spike_data_',num2str(trial)];
        spike_data = EZ.(spikeDataName);
        extractionHeaderName = ['extraction_header_',num2str(trial)];
        extraction_header = EZ.(extractionHeaderName);
        extractionDataName = ['extraction_data_',num2str(trial)];
        extraction_data = EZ.(extractionDataName);
        trialStatus = trial_header.success;
        
        %Get Target data
        targetID = trial_header.target_id;
        targetLoc = 1000*trial_header.target.distance*[cosd((trial_header.target_id-1)*45), sind((trial_header.target_id-1)*45)];
        
        %Get decoded position/vel
        pos = extraction_data.pos.*1000;
        vel = extraction_data.vel.*1000;
        send_time_vel = extraction_header.send_time_vel;
        send_time_pos = extraction_header.send_time_pos;
            
            trialStartTimeComputer = trial_header.task_state.time(1);
            posTime = round(1000*(send_time_pos-trialStartTimeComputer));
            velTime = round(1000*(send_time_vel-trialStartTimeComputer));
                        
       %Check state times
        stateNames = {'SessionStart','Center','HoldA','React','Move','Hold','InterTrial','TrialEnd'};
        trialStateNames = trial_header.beh_event.type;
        numTrialStates = numel(trialStateNames);
        stateTransitions = zeros(2,numTrialStates);
        for trialState = 1:numTrialStates
           stateTransitions(1,trialState) = find(cellfun(@(x) strcmp(trialStateNames{1,trialState},x), stateNames)); 
        end
        %Align to state transitions to start, convert to ms, round
        stateTransitions(2,:) = round(1000*(trial_header.beh_event.time-trial_header.beh_event.time(1,1)))';

        

        reactTime = stateTransitions(2,find(stateTransitions(1,:)==4));
        moveTime = stateTransitions(2,find(stateTransitions(1,:)==5));
        holdTime = stateTransitions(2,find(stateTransitions(1,:)==6));
       
        %Plot decoded pos/vel with state transtions
        figure
            subplot(4,1,1); hold on
                plot(posTime,pos(1,:))
                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('x')
            subplot(4,1,2); hold on
                plot(posTime,pos(2,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('y')
            subplot(4,1,3); hold on
                plot(velTime,vel(1,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('Vx')
            subplot(4,1,4); hold on
                plot(velTime,vel(2,:))
                                ax = gca;
                line([reactTime,reactTime],[ax.YLim])
                line([moveTime,moveTime],[ax.YLim])
                line([holdTime,holdTime],[ax.YLim])
                ylabel('Vy')
        
        