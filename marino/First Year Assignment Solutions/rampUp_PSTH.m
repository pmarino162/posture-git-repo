% Assignment: Peristimulus Time Histogram (PSTH)
%% Setup
close all, clearvars -except Data
if ~exist('Data','var')||~isfield(Data(1).TrialData,'timeMoveOnset')
    load('Ike20120905_processed_original.mat');
end
o = [Data.Overview];
successData = Data([o.trialStatus] == 1);
successData(end) = []; %remove last trial (no spikes present)

startT = tic;
dt = 20;
binT = -500:dt:500;
fr = NaN(length(successData),...
    length(successData(1).TrialData.spikes),length(binT)-1);
for trl = 1:length(successData)
    t = successData(trl).TrialData.timeMoveOnset+binT;
    for nrn = 1:length(successData(trl).TrialData.spikes)
        spkT = successData(trl).TrialData.spikes(nrn).timestamps;
        for bn = 1:(length(binT)-2)
            fr(trl,nrn,bn) = sum(spkT>=t(bn)&spkT<t(bn+1));
        end
        fr(trl,nrn,end) = sum(spkT>=t(bn)&spkT<=t(bn+1));
    end
end
endT = toc(startT(1));

%% Part 1 (Trial Averaged PSTH)
nrn = 10;
meanFR = squeeze(mean(fr,1));
figure, bar(binT(1:(end-1))+10,meanFR(nrn,:)/(dt/1000))
title(sprintf('Neuron %i',nrn),'fontsize',18)
xlabel('Time from Movement Onset (ms)','fontsize',14)
ylabel('Firing Rate (Hz)','fontsize',18)

%% Part 2 (Target Averaged PSTH)
nrn = 172;
td = [successData.TrialData];
ang = [td.reachAngle];
[uniqAng,~,angIdx] = unique(ang);
figure
plotPos = [6 3 2 1 4 7 8 9];
ax = NaN(length(uniqAng),4);
for i = 1:length(uniqAng)
    subplot(3,3,plotPos(i))
    bar(binT(1:(end-1))+10,squeeze(mean(fr(angIdx==i,nrn,:),1)/(dt/1000)))
    ax(i,:) = axis(gca);
    switch plotPos(i)
        case 2; title(sprintf('Neuron %i',nrn),'fontsize',18);
        case 4; ylabel('Firing Rate (Hz)','fontsize',14);
        case 8; xlabel('Time from Movement Onset (ms)','fontsize',14);
    end
    set(gca,'fontsize',12)
end
[~,idx] = max(ax(:,end));
for i = 1:length(uniqAng); subplot(3,3,plotPos(i)), axis(ax(idx,:)); end