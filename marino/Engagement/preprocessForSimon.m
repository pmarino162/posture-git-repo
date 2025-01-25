clear; clc; close all;

%% Load coloarmap
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
 
for date = {'20200610','20200611','20200612','20200615'}
    %% Load Data
    if strcmpi(date,'20200610')
        load('D:\Animals\Earl\2020\06\20200610\00_delayedCenterOut\Earl20200610_00_delayedCenterOut_translated.mat')
    elseif strcmpi(date,'20200611')
        load('D:\Animals\Earl\2020\06\20200611\00_delayedCenterOut\Earl20200611_00_delayedCenterOut_translated.mat')
    elseif strcmpi(date,'20200612')
        load('D:\Animals\Earl\2020\06\20200612\00_delayedCenterOut\Earl20200612_00_delayedCenterOut_translated.mat')
    elseif strcmpi(date,'20200615')
        load('D:\Animals\Earl\2020\06\20200615\00_delayedCenterOut\Earl20200615_00_delayedCenterOut_translated.mat')
    end
    
    %% Preprocess Data
    getSpikes = false;
        getSorts = false;
        exclCh =  [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = true;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = false;
        centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = false;
        forceSetup = 'EL';
    getAlg = true;
    getEye = true;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = false;
    checkDecode = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];

    Data = getDataStruct20211210(Data,'getSpikes',getSpikes,'getSorts',getSorts,'exclCh',exclCh,'exclZero',exclZero,...
        'getMarker',getMarker,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'removeBadPhasespace',removeBadPhasespace,'centerMarker',centerMarker,...
        'getForce',getForce,'forceSetup',forceSetup,'getAlg',getAlg,'getEye',getEye,'getWaveforms',getWaveforms,'getLFP',getLFP,'getBB',getBB,'getKin',getKin,...
        'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'inclStateTable',inclStateTable,'exlcTrials',exclTrials);  
    
% %% Make sure it worked
%     Data = Data([Data.trialStatus]==1);
%     rmTrials = [];
%     for trial = 1:size(Data,2)
%         if isempty(Data(trial).eye.position)
%             rmTrials = [rmTrials,trial];
%         end
%     end
%     Data(rmTrials) = [];
%     
%     binWidth = 25;
%     trajFields = {'marker','eyePos','pupil'};
%     condFields = {{'target','targetData','targetID'}};
%     
%     trialInclStates = struct('trialName','','inclStates',[]);
%         trialInclStates(1).trialName = {'HC_CenterOut'};
%         trialInclStates(1).inclStates = {{'state','Step 2','first',0},{'state','Success with Reward','first',0}};
% 
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

% %% Plot 
%     figure
%     for target = 1:8
%         traj = trajStruct([trajStruct.target]==target).avgMarker.traj;
%         plot(traj(:,1),traj(:,2),'Color',tcmap(target,:))
%         hold on
%     end
%     
%     figure
%     for target = 1:8
%         traj = trajStruct([trajStruct.target]==target).avgEyePos.traj;
%         plot(traj(:,1),traj(:,2),'Color',tcmap(target,:))
%         hold on
%     end
    
%     %% Plot pupil
%     figure
%     plotTime = 0;
%     for trial = 1:size(Data,2)
%         time = Data(trial).eye.time;
%         time = time + plotTime;
%         plotTime = time(end);
%         pupil = Data(trial).eye.pupil(:,1);
%         if ~isempty(pupil)
%             plot(time,pupil)
%             hold on
%         end
%     end
    
%% Save Data
    filepath = 'D:\Animals\Earl\2020\06\PreprocessedForSimon\'
    filename = [filepath,'Earl',date{1,1},'.mat']
    save(filename,'Data')

end