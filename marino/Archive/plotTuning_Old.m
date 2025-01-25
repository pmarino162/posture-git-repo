function plotTuning(Data)

%% Specify Time Window
    timeWindow = [720 900];
    windowLength = timeWindow(2)-timeWindow(1);
    
%% Get Only Successful Trials
    Data = Data([Data.trialStatus]==1);
    
%% Get Tuning Data
    numChannels = size(Data(1).spikes,2);
    targetData = [Data.targetData];
    targetID = [targetData.targetID];
    numTargets = length(unique(targetID));
    
   for target = 1:numTargets
      tempData = Data(targetID==target);
      for trial = 1:size(tempData,2)
          binTimes = tempData(trial).spikes(1).binTimes;
          binMask = binTimes >= timeWindow(1) & binTimes <= timeWindow(2);
          for channel = 1:numChannels
            trialFR = (sum(tempData(trial).spikes(channel).allSorts(binMask))./windowLength).*1000;
            tuningData.channels(channel).allData(target).allTargetData(trial) = trialFR;
          end
      end
      
      for channel = 1:numChannels
        allTargetData = tuningData.channels(channel).allData(target).allTargetData;
        meanChannelData = mean(allTargetData);
        stdChannelData = std(allTargetData);
        numObs(channel,target) = size(allTargetData,2);
        tuningData.channels(channel).Means(target) = meanChannelData;
        tuningData.channels(channel).SD(target) = stdChannelData; 
      end
   end
   
   %Fit Tuning Curves
    targetAngles = transpose(0:45:315);
    x = [ones(8,1),sind(targetAngles),cosd(targetAngles)];

    for channel = 1:numChannels
     minNumObs = min(numObs(channel,:));
        y = []; x = [];
        for target = 1:8
            allChData = tuningData.channels(channel).allData(target).allTargetData(1:minNumObs)';
            targData = [ones(minNumObs,1),ones(minNumObs,1)*sind(targetAngles(target)),ones(minNumObs,1)*cosd(targetAngles(target))];
            y = vertcat(y,allChData);
            x = vertcat(x,targData);
        end
%         y = tuningData(1).channels.Means(channel,:)';
        [B,bint,r,rint,stats] = regress(y,x);
%             B = x\y; 
        b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); p = stats(3);
        tuningData.channels(channel).MD = sqrt(b1.^2 + b2.^2);
        tuningData.channels(channel).PD= atan2d(b1,b2);
%         if atan2d(b1,b2) < 0
%             tuningData.channels(channel).PD= 360+atan2d(b1,b2);
%         end
        tuningData.channels(channel).b0 = b0;
        tuningData.channels(channel).p = p;
    end
    
    
%% Plot Data for each Time Bin 
    f=figure('Position',[0 0 600 600])
    f.Name = [date,'Vector']
    for channel = 1:numChannels
        PD = tuningData.channels(channel).PD;
        modDepth = tuningData.channels(channel).MD;
%         if ismember(channel,inclCh)
%             quiver(0,0,modDepth*cosd(PD),modDepth*sind(PD),'LineWidth',2,'Color','r')
%         else
            quiver(0,0,modDepth*cosd(PD),modDepth*sind(PD),'LineWidth',2,'Color','b')
%         end
        hold on
    end
    sgtitle([date,' Channels'])
%     saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-5-19')
%     xlim([-3 3])
%     ylim([-3 3])
% 
    f=figure('Position',[0 0 1200 800])
    f.Name = [date,'Tuning']
    for channel = 1:numChannels   
        subplot(10,10,channel)
        avgFR = tuningData.channels(channel).Means;
        stdFR = tuningData.channels(channel).SD;
        PD = tuningData.channels(channel).PD;
        modDepth = tuningData.channels(channel).MD;
        b0 = tuningData.channels(channel).b0;
        p = tuningData.channels(channel).p;
        Bfit = [b0;modDepth];
        x = [ones(8,1),cosd(targetAngles-PD)];
        cosFit = x*Bfit;
        if p < 0.05
%             if ismember(channel,inclCh)
%                 plot(targetAngles,cosFit,'Color','r')
%             else
                plot(targetAngles,cosFit,'Color','b')
%             end
        else
%             if ismember(channel,inclCh)
%                 plot(targetAngles,cosFit,'--','Color','r')
%             else
                plot(targetAngles,cosFit,'--','Color','b') 
%             end
        end
        hold on
        errorbar(targetAngles,avgFR,stdFR,'.k')
        plot(targetAngles,avgFR,'.','MarkerSize',2);
        title(['Ch ',num2str(channel)])
    end
    sgtitle([date,' Channels'])
%     saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-5-19')


end