% GPFA decoder online calibration script

%% Setup
setupPath({'rigAnalysisCode','degenhart','projects/EnergyLandscape'})

%% Set parameters
% Parameters that should be changed daily
exclCh = [13 31 10];
workspaceCenter = [105 335];

[P] = initGPFADecoderParams('exclCh',exclCh,'workspaceCenter',workspaceCenter);

% Define perpendicular scale values
perpScaleVal = [0.001 0.25 0.5 1 1];

% Define file paths
% baseDir = 'E:\Animals';
subject = 'Earl';
% dn = now;
% yr = datestr(dn,'yyyy');
% mo = datestr(dn,'mm');
% dataset = datestr(dn,'yyyymmdd');
% dsPath = fullfile(baseDir,subject,yr,mo,dataset);
dsPath = 'E:\Animals\Earl\2020\01\20200124';
dataset = '20200124';
calibFilenames = { ...
    '02_observation_02', ...
    '03_gradualTraining', ...
    '03_gradualTraining', ...
    '03_gradualTraining', ...
    '03_gradualTraining'};

%% Create trialsToUse Struct
trialsToUse = struct('decoderID',[],'calibFilename',[],'allTrials',[],'artifactTrials',[],'forceBarFailures',[],'goodTrials',[],'targetBalancedTrials',[],'trialsToUse',[]);

%% Subselect Trials To Remove Artifact 
% Setup paths
pN = fullfile(dsPath,calibFilenames{P.decoderID + 1});
fN = [subject dataset '_' calibFilenames{P.decoderID + 1} '_SI_translated.mat'];
%Subselect Trials
[trialsToUse] = subselectToRemoveArtifact(pN,fN,P.decoderID + 1,calibFilenames{P.decoderID + 1},trialsToUse)
trialsToUseOut = trialsToUse(P.decoderID + 1).trialsToUse;

%% Calibrate intuitive decoder -- Run x5

% Define calibration parameters
P.removeJaggedDim = false;
P.jaggedDim = [];
P.binwidth = 45;
P.nBins = 7;
successOnly = false;

tauBinThresh = 0;
kinOption = 'endOfTrial'; % 'maxDist', 'endOfTrial', 'targetRamp'
scaleFactor = 1;
excludeLatents = [];
xSpec = 'xorth';

% Setup paths
pN = fullfile(dsPath,calibFilenames{P.decoderID + 1});
fN = [subject dataset '_' calibFilenames{P.decoderID + 1} '.mat'];

% Run calibration
[P,F] = calibrateGPFADecoder(P,'kinOption',kinOption, ...
    'tauBinThresh',tauBinThresh,'tauScaleFactor',scaleFactor,...
    'excludeLatents',excludeLatents,'xSpec',xSpec, ...
    'pathName',pN,'fileName',fN,'successOnly',successOnly,'trialsToUseOut',trialsToUseOut);
P.perpendicularScale = perpScaleVal(P.decoderID);

% Save decoder parameters and figures
saveFigs = true;

saveGPFADecodingParams(P)

if saveFigs
    saveDir = fullfile(P.analysisDir,P.decoderName);
    mkdir(saveDir)
    saveFigurePDF(F,saveDir)
end
close all

%% Run GPFA

% Define file paths
subject = 'Earl';
dataDir = '13_condGridTask_04';%'04_condGridTask_01';
gpfaRunID = 13;

% Setup file paths
baseDir = 'D:\Data';
dn = now;
yr = datestr(dn,'yyyy');
mo = datestr(dn,'mm');
dataset = datestr(dn,'yyyymmdd');
dsPath = fullfile(baseDir,subject,yr,mo,dataset);
dataDir = fullfile(dsPath,dataDir);

% Define parameters
xSpec = 'xorth';
binWidth = 45;
causal = true;
causalMode = 'rt';


[G,E,F1] = runGridTaskGPFA(subject,dataset, ...
    'exclCh',exclCh,'xSpec',xSpec,'binWidth',binWidth, ...
    'causal',causal,'causalMode',causalMode,'gpfaRunID',gpfaRunID, ...
    'dataDir',dataDir);

opt_nAvg = 1;
objectiveFn = 'standard';

[E.OptProj,F2] = identifyOptimizedProjection(G,E,'selectProj',false, ...
    'opt_nAvg',opt_nAvg,'plotObjFn',true,'objectiveFn',objectiveFn);

%saveFigurePDF([F1 F2],saveDir)
%close([F1 F2])

%% Generate decoder for rotated mapping
projID = 1;
E.invSec = true;
nAvg = 1;
alignMode = 'start';

[P,F3,Gtrans] = fitRotatedGPFAMapping(G,E,projID,'nAvg',nAvg, ...
    'alignMode',alignMode);

%% Plot the GPFA Cloud
% 
% [F4,cloudInfo] = plotRotatedNeuralCloud([Gtrans Gtrans Gtrans],E.OptProj,'saveFigs',0,...
%     'saveName',E.saveNameBase);
[F4,cloudInfo] = plotRotatedNeuralCloud([Gtrans],E.OptProj(projID),'saveFigs',0,...
    'saveName',E.saveNameBase);



% Save decoder parameters and figures
saveFigs = true;
%P.decoderName(end-3:end) = 'T1T5';
P.decoderName(end-7:end) = 'T1T5_inv';
saveGPFADecodingParams(P)

% Save the figures about the optimization and the latents
saveDir = fullfile(E.analysisPath,[P.decoderName '_gpfaRunID_' num2str(gpfaRunID)]);
mkdir(saveDir)

if saveFigs
    %F = [F1 F2 F3 F4];
    F = [F3 F4];
    saveDir = fullfile(E.analysisPath,[P.decoderName '_gpfaRunID_' num2str(gpfaRunID)]);
    mkdir(saveDir)
    saveFigurePDF(F,saveDir)
end
%close F3 F4
saveFigurePDF([F1 F2],saveDir)
close all

%% Find the intermediate target locations (based on the two-target rotated mapping data)

[sF,intTargPos] = findIntermediateTargetLocations