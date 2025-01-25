function [kinStruct,neuronStruct,Data] = generateNeuralAndKinData(posLag,velLag)

    %% Set parameters 
    numCondTrials = 20;
    numNeurons = 96;
    ts = 25; %ms
    postureList = [1:3];
    postureLocations = [-50,-50; 50,-50; 50,50; -50,50];
%     postureList = 1:7;
%     postureLocations = [-120,40; -40,40; 40,40; 120,40; 0,0; -40,-40; 40,-40]; 
    targetList = 1:8;
    targetDist = 100; %mm
    reachTime = 700; %ms
    holdTime = 300; %ms
    
    profileMode = 'bellVelStatPos';
    

    %% Create 1D position and velocity profiles
    time = 0:ts:reachTime+holdTime;
    numReachSteps = sum(time<reachTime);
    numHoldSteps = sum(time>=reachTime);
    numSteps = size(time,2);
    
    switch profileMode
        case 'constVel'
            %Constant velocity
            velProf = [(targetDist/reachTime)*ones(1,numReachSteps),zeros(1,numHoldSteps)]; %mm/ms
            posProf = zeros(1,numSteps);
            for i = 2:numSteps
               posProf(i) = posProf(i-1) + velProf(i-1)*ts; 
            end
        case 'bellVel'
            %Bell-shaped velocity profile
            peakSpeed = 0.5; %mm/ms
            peakSpeedTime = 300;%ms
            a = pi/((targetDist/peakSpeed)^2);
            velProf = [peakSpeed*exp(-a*(time(time<reachTime)-peakSpeedTime).^2),zeros(1,numHoldSteps)]; %mm/ms
            posProf = zeros(1,numSteps);
            for i = 2:numSteps
               posProf(i) = posProf(i-1) + velProf(i-1)*ts; 
            end
        case 'constVelStatPos'
            %Constant velocity; static position
            velProf = [(targetDist/reachTime)*ones(1,numReachSteps),zeros(1,numHoldSteps)]; %mm/ms
            posProf = zeros(1,numSteps);
        case 'bellVelStatPos'
            %Bell-shaped velocity profile; static position
            peakSpeed = 0.5; %mm/ms
            peakSpeedTime = 300;%ms
            a = pi/((targetDist/peakSpeed)^2);
            velProf = [peakSpeed*exp(-a*(time(time<reachTime)-peakSpeedTime).^2),zeros(1,numHoldSteps)]; %mm/ms
            posProf = zeros(1,numSteps);
    end
    
    figure
    subplot(2,1,1)
        plot(time,velProf)
        xlabel('time (ms)')
        ylabel('speed (mm/ms)')
    subplot(2,1,2)
        plot(time,posProf)
        xlabel('time (ms)')
        ylabel('distance (mm)')
     close
        
%% Get position and velocity trajectories
    kinStruct = struct('posture',[],'target',[],'time',[],'position',[],'velocity',[]);
    structInd = 1;
    for posture = postureList
        for target = targetList
            theta = (target-1)*45;
            position = zeros(numSteps,2);
            velocity = zeros(numSteps,2);
            startPos = postureLocations(posture,:);   
            for i = 1:numSteps
                position(i,:) = startPos + [posProf(i)*cosd(theta), posProf(i)*sind(theta)];
                velocity(i,:) = [velProf(i)*cosd(theta), velProf(i)*sind(theta)];
            end
            kinStruct(structInd).posture = posture;
            kinStruct(structInd).target = target;
            kinStruct(structInd).time = time;
            kinStruct(structInd).position = position;
            kinStruct(structInd).velocity = velocity;
            structInd = structInd + 1;
        end
    end
    
%     for posture = postureList
%         for target = targetList
%             time = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).time;
%             position = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).position;
%             velocity = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).velocity;
%             
%             figure
%             subplot(2,2,1)
%                 plot(time,position(:,1)); xlabel('time (ms)'); ylabel('x (mm)')
%             subplot(2,2,2)
%                 plot(time,position(:,2)); xlabel('time (ms)'); ylabel('y (mm)')
%             subplot(2,2,3)
%                 plot(time,velocity(:,1)); xlabel('time (ms)'); ylabel('V_x (mm/ms)')
%             subplot(2,2,4)
%                 plot(time,velocity(:,2)); xlabel('time (ms)'); ylabel('V_y (mm/ms)')
%             sgtitle([profileMode,' T',num2str(target),' P',num2str(posture)])
% %             saveas(gcf,[dirStr,'Kinematics\',profileMode,' T',num2str(target),' P',num2str(posture),' kinematics.jpg'])
%             close 
%         end
%     end
    
%% Generate Neuron Properties
    neuronStruct = struct('bn',[],'Bv',[],'Bp',[],'b0',[]);
%     meanbn = 40*(ts/1000);
%     stdbn = 1*(ts/1000);
%     meanMagBv = 40*(ts/1000);
%     stdMagBv = meanMagBv/5;
%     meanMagBp = .2*(ts/1000);
%     stdMagBp = meanMagBp/5;
%     meanb0 = 22.5*(ts/1000);
%     stdb0 = 5*(ts/1000);
    meanbn = 25*(ts/1000);
    stdbn = meanbn/5;
    meanMagBv = 45*(ts/1000);
    stdMagBv = meanMagBv/5;
    meanMagBp = .08*(ts/1000);
    stdMagBp = meanMagBp/5;
    meanb0 =30*(ts/1000);
    stdb0 = 2*(ts/1000);
    for neuron = 1:numNeurons
       %Speed coefficient 
       bn = normrnd(meanbn,stdbn); 
       %Preferred vel direction
       PD = randn(1,2);
       PD = PD./norm(PD);
       %Preferred pos direction
       PP = randn(1,2);
       PP = PP./norm(PP);
       %Vel mod depth
       Mv = normrnd(meanMagBv,stdMagBv);
       %Pos mod depth
       Mp = normrnd(meanMagBp,stdMagBp);
       %Baseline FR
       b0 = normrnd(meanb0,stdb0); 
       %Model params
       neuronStruct(neuron).bn = bn;
       neuronStruct(neuron).Bv = Mv.*PD;
       neuronStruct(neuron).Bp = Mp.*PP;
       neuronStruct(neuron).b0 = b0;
    end

% %% PDs
%     figure('Position',[0 0 600 600])
%     color = [1 0 0];
%     for neuron = 1:numNeurons
%         Bv = neuronStruct(neuron).Bv;
%         quiver(0,0,Bv(1),Bv(2),'LineWidth',2,'Color',color)
%         hold on
%     end
%     sgtitle('Pref. Dirs scaled by Mod Depth')

%% Generate Data Struct 
    stateData = struct('stateNames',{},'stateTransitions',[]);
    marker = struct('time',[],'position',[],'velocity',[]);
    Data = struct('trialNum',[],'trialName','','trialStatus',[],'targetData',[],'conditionData',[],'stateData',stateData,'marker',marker);
    sortList = struct('channel',[],'sort',[]);
    for i = 1:numNeurons
       sortList(i).channel = i;
       sortList(i).sort = 'all';
    end
    posLagBins = round(posLag/ts);
    velLagBins = round(velLag/ts);
    structInd = 1;
    trialNum = 0;
    for posture = postureList
        for target = targetList
            for trial = 1:numCondTrials
                %Create Marker Data
                time = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).time;
                velocity = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).velocity;
                position = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).position; 
                %Create Neural Data
                allChannelSpikeBins = zeros(size(time,2),numNeurons);
                for neuron = 1:numNeurons
                    bn = neuronStruct(neuron).bn;
                    Bv = neuronStruct(neuron).Bv;
                    Bp = neuronStruct(neuron).Bp;
                    b0 = neuronStruct(neuron).b0;
%                     noise = normrnd(0,.1,numSteps,1);
                    velNorm = vecnorm(velocity,2,2);
%                     allChannelSpikeBins(:,neuron) = bn*velNorm + velocity*Bv' + position*Bp' + b0 + noise;
                    expectedSpikes = NaN(size(time,2),1);
                    for i = 1:size(expectedSpikes,1)-max(velLagBins,posLagBins)
                        expectedSpikes(i,1) = bn*velNorm(i+velLagBins,1) + velocity(i+velLagBins,:)*Bv' + position(i+posLagBins,:)*Bp' + b0;
                    end
                    expectedSpikes(expectedSpikes<0) = 0;
                    if any(expectedSpikes<0)
                        fprintf('exp Spikes < 0 \n')
                    end
                    poissonSpikes = poissrnd(expectedSpikes);
                    
                    allChannelSpikeBins(:,neuron) = poissonSpikes;
                end
                %Add to Data Struct
                Data(structInd).trialNum = trialNum + 1;
                Data(structInd).trialName = 'Center Out';
                Data(structInd).conditionData.postureID = posture;
                Data(structInd).targetData.targetID = target;
                Data(structInd).stateData(1).stateNames = {'Reach','Hold','Success'};
                Data(structInd).stateData(1).stateTransitions = [1,2,3;0,reachTime,reachTime+holdTime];
                Data(structInd).spikes.binTimes = time;
                Data(structInd).spikes.allChannelSpikeBins = allChannelSpikeBins;
                Data(structInd).spikes.sortList = sortList;
                Data(structInd).marker.time = time;
                Data(structInd).marker.position = position;
                Data(structInd).marker.velocity = velocity;
                structInd = structInd + 1;
            end
        end
    end
    



end