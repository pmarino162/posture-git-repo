%Problem 1
figure;
time = table2array(data(:,1));
bb = table2array(data(:,2));

plot(time,bb,'b-','LineWidth',2);
xlabel('time (s)')
ylabel('bb (uV)')

%Problem 3
spikesBool = logical(table2array(data(:,3)));
spikeTimes = time(spikesBool,1);
bbAtSpikes = bb(spikesBool,1);

figure; hold on;
plot(time,bb,'b-','LineWidth',2);
plot(spikeTimes,bbAtSpikes,'r.','MarkerSize',15)
xlabel('time (s)')
ylabel('bb (uV)')

%Problem 4
duration = time(end) - time(1);
numSpikes = size(spikeTimes,1);
spikeRate = numSpikes/duration;

%Problem 5
figure; hold on;
for spike = 1:numSpikes
    %Get spike time
    curSpikeTime = spikeTimes(spike,1);
    %Get time range
    tStart = curSpikeTime - 2/1000;
    tEnd = curSpikeTime + 1/1000;
    %Select bb snippet
    timeMask = time > tStart & time < tEnd;
    curWaveform = bb(timeMask,1);
    plot(curWaveform);
end
xlabel('sample time')
ylabel('bb (uV)')

%Problem 6

%Filter bb
% d1 = designfilt("lowpassiir",FilterOrder=12, ...
%     HalfPowerFrequency=0.15,DesignMethod="butter");

fc = 300;
fs = 1.958*10^4;
[b,a] = butter(4,fc/(fs/2));
filtBB = filtfilt(b,a,bb);
%Plot filtBB
figure; hold on;
plot(time,filtBB);
%Threshold filtBB
deBiasedBB = bb-filtBB;
figure; hold on
plot(time,deBiasedBB,'b');
newSpikeTimes = time(deBiasedBB > 50);


plot(newSpikeTimes,deBiasedBB(deBiasedBB > 50),'r.','MarkerSize',15)
