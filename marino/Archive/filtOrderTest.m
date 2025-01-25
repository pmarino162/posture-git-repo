clear;
clc;
clf;

%Filter test

%% Load data 
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
%% Preprocess data
    getMarker = false;
    getKin = false;
    checkPhasespaceSync = false;
    getForce = false;
    getAlg = false;
    getEye = false;
    getSpikes = false;
    getSorts = false;
    getSmoothedFR = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    inclStateTable = false;
    removeBadPhasespace = true;
    binWidth = 1;
    exclCh = [];
    exclTrials = [];
    exclZero = true; %Exclude channels that only contain sort 0
    centerMarker = true; %Subracts workspace center from marker pose if true
    forceSetup = '';
%     assignopts(who,varargin);
    procData = preprocessDataMarino(Data,'CHECKPHASESPACESYNC',checkPhasespaceSync,'DROPPEDMARKERTRIALS','all','CHECKDECODE',false);

%% Get raw data for a trial 
    trial = 1;
    rawPositions = procData(trial).TrialData.Marker.rawPositions;
    time = rawPositions(:,6)';
    position = rawPositions(:,2:3);
    %Flip Y
    position(:,2) = -1.*position(:,2);

%% Design filter
    time = rawPositions(:,6)';
    fs =  1/(mode(diff(time))/1000);                       %sampling frequency (Hz)
    Wn = 15/fs;                                            %cutoff frequency (normalized)
    [b,a] = butter(4,Wn);  
    
%% Filter position; get velocity
    position1 = filtfilt(b,a,position);
    %Interpolate up to 1000Hz
    interpTime = ceil(time(1)):1:floor(time(end));
    position3 = interp1(time,position1,interpTime);
    %Get velocity, interpolate to position times
    dp = diff(position3);
    dt = diff(interpTime);
    velocity1 = dp./dt';
    velocityTime1 = interpTime(1:end-1)+dt/2;
%     velocity1 = interp1(velocityTime1,velocity1,interpTime);

    
    position2 = filtfilt(b,a,position);
    dp = diff(position2);
    dt = diff(time);
    velocity2 = dp./dt';
    velocityTime2 = time(1:end-1)+dt/2;
    velocity2 = interp1(velocityTime2,velocity2,interpTime);
    position2 = interp1(time,position2,interpTime);
%Compare
figure
plot(interpTime,position1)
hold on
plot(interpTime,position2)

figure
plot(velocityTime1,velocity1)
hold on
plot(interpTime,velocity2)