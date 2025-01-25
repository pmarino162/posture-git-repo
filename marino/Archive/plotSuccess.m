function plotSuccess(Data)


windowSize = 50;
trialNum = [Data.trialNum];
numTrials = length(trialNum);
trialStatus = [Data.trialStatus];

Decoder = [Data.Decoder];
decoderName = vertcat(Decoder.name);

decoderChangeLoc = [];

for trial = 2:numTrials
    if ~strcmpi(decoderName(trial,:),decoderName(trial-1,:))
        decoderChangeLoc = [decoderChangeLoc,trial];
    end
end

movAvg = (movmean(trialStatus,windowSize)).*100;

plot(trialNum,movAvg);
hold on
for i = 1:length(decoderChangeLoc)
    xline(decoderChangeLoc(i))
end
xlabel('Trial')
ylabel('Success Rate (%)')
ylim([0 100])

end