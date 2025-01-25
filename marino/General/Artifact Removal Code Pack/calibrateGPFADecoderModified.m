function [P,F] = calibrateGPFADecoder(P,varargin)
% Calibrate GPFA decoder
%
% [P] = calibrateGPFADecoder(P)

%--------------------------------------------------------------------------
% Parse optional arguments

% Optional arguments
removeOldGPFAResults = false;
tWin = 100;         % Time window (in ms) used during calibration
tauBinThresh = 2;   % Threshold for excluding short timescale latents
fileName = [];      % Data filename
pathName = [];      % Path to data
kinOption = 'endOfTrial'; % Method to generate kinematic training data
nStartSamp = 3;  % Number of start samples to use
nTargSamp = 5;  % Number of target/movement samples to use
kinScaleFactor = 1;
tauScaleFactor = 1;
plotPrecomp = false;
excludeLatents = [];
xSpec = 'xorth';
analysisMode = 'online';
successOnly = true; % Specify whether or not to use successful trials only when calibrating
analysisDir = 'artifactAnalysis';
trialsToUseOut = [];
assignopts(who,varargin);

%--------------------------------------------------------------------------
% Setup and load data

% Prompt user to select data to calibrate
if isempty(fileName) && isempty(pathName)
    [fileName,pathName] = uigetfile;
end
f = fileName;
p = pathName;
dataPath = fullfile(p,f);

% Load data, preprocess, and convert to trajectory data object
% If running in online mode, use UDP spikes.  If running in offline mode, 
% use TDT spikes.
switch analysisMode
    case 'online'
        optArgs = {'CHECKPHASESPACESYNC',false,'GETUDP',true};
    case 'offline'
        optArgs = {'CHECKPHASESPACESYNC',false};
end
Data = load(dataPath);
Data = Data.Data;
procData = preprocessData(Data,optArgs{:});
if ~isempty(trialsToUseOut)
procData = procData(trialsToUseOut);
end
T = convertToClassFormat(procData);

%validStates = {'Start','Touch Bar Reach','Touch Bar Hold',

TD = getTrajectoryDataModified(T);
TD = TD.normalize;

% Set up analysis and save paths
P.subject = TD(1).subject;
P.date = datestr(TD(1).date,'yyyymmdd');
%pathParts = split(p,filesep);
sepInd = find(ismember(p,filesep));
pathOnset = sepInd(end) + 1;
pathOffset = length(p);
P.dataDir = p(pathOnset:pathOffset);
P.baseDir = fullfile(p(1:sepInd(end)));
P.saveNameBase = [P.subject P.date '_' P.dataDir];
P.analysisDir = fullfile(P.baseDir,analysisDir);

% Create analysis directory (if needed)
mkdir(P.analysisDir)

%--------------------------------------------------------------------------
% Run GPFA

% Remove old GPFA results.  This is done for debugging purposes so that
% previously-computed parameters aren't re-used.
if removeOldGPFAResults
    [~,~,~] = rmdir(fullfile(P.analysisDir,'mat_results'),'s');
end

% Increment decoder ID.  This is also the ID used to save GPFA results.
P.decoderID = P.decoderID + 1;
P.decoderType = 'GPFA';
P.decoderName = sprintf('MkHstDec_GPFA_%s%s_Intuitive_%0.2d', ...
    P.subject,P.date,P.decoderID);

% Set parameters for causal GPFA based on bin width if not specified
if P.causalGPFA && strcmp(P.causalGPFAMode,'rt') && isempty(P.nBins)
    P.Tmax = round(1000/P.binwidth);
    P.nBins = round(300/P.binwidth);
end

% Run GPFA
G = gpfaAnalysis(TD,P.exclCh,P.decoderID,P.analysisDir, ...
    'validSort',P.validSort,'binWidth',P.binwidth,'onsetMode',P.onsetMode, ...
    'removeJaggedDim',P.removeJaggedDim,'jaggedDim',P.jaggedDim, ...
    'saveResults',false,'causal',P.causalGPFA,'causalMode',P.causalGPFAMode, ...
    'segLength',P.segLength,'causal_Tmax',P.Tmax,'causal_nBins',P.nBins, ...
    'tauScaleFactor',tauScaleFactor,'plotPrecomp',plotPrecomp);

% Plot GPFA results -- Note, might need to remove the 'ColMat' part of the
% following functions if using old MATLAB versions (or if the repository
% isn't up to date)
% C = el.defineTaskColormap('centerOut');
% F(1) = plotLatentDimVsTime(G, ...
%     'ColMat', C, ...
%     'saveNameBase',P.saveNameBase);
% F(2) = plotLatentDimVsTime(G, ...
%     'xspec','xorth', ...
%     'ColMat', C, ...
%     'saveNameBase',P.saveNameBase);

%--------------------------------------------------------------------------
% Get calibration data -- use distance from baseline method

fprintf('\nGetting calibration data ... ')

% Get successful trials ONLY (NOTE: this could be changed to use all
% attempted trials.
if successOnly
    successMask = logical([G.TD.successful]);
    G.TD = G.TD(successMask);
end
nTrials = length(G.TD);

% Loop over trials, get data at the start to define baseline
nTs = floor(tWin/P.binwidth);
P.tWin = tWin;
P.nTs = nTs;
C = repmat(struct(),nTrials,1);
for i = 1:nTrials
    % Get start (center) and end (target) positions
    startPos = G.TD(i).startPos;
    endPos = G.TD(i).targPos;
    
    % Find task timing.  Get this from the 'stateOnset' field and divide by
    % the movement onset to find the nearest bin.
    onset = round(G.TD(i).stateOnset/P.binwidth) + 1;
    
    % Get GPFA data and truncate to match trajectory
    y = G.TD(i).GPFA.(xSpec);
    y = y(:,1:end);
    nSamp = size(y,2);
    
    % Get binned spike counts
    z = G.TD(i).GPFA.y;
    z = z(:,1:end);

    % Average baseline activity
    yB = mean(y(:,1:nTs),2);
    zB = mean(z(:,1:nTs),2);
    mask_base = false(1,nSamp);
    mask_base(1:nTs) = true;
    
    % Add to training data structure
    C(i).trialID = G.TD(i).trialID;
    C(i).successful = logical(G.TD(i).successful);
    C(i).nSamp = nSamp;
    C(i).startPos = startPos;
    C(i).endPos = endPos;
    C(i).onset = onset;
    C(i).yB = yB;   % baseline latent state
    C(i).zB = zB;   % baseline spike count
    C(i).y = y;     % latent state
    C(i).z = z;     % spike counts
    C(i).d_z = [];  % distance from baseline (latent state)
    C(i).d_y = [];  % distance from baseline (spike count)
    C(i).mask_base = mask_base;  % baseline sample mask
    C(i).mask_move = [];         % movement sample mask
    C(i).max_ind_z = [];         % index of max distance (spike count)
    C(i).max_ind_y = [];         % index of max distance (latent state)
    C(i).Ystart = [];
    C(i).Yend = [];
    C(i).Y = [];
    C(i).Xstart = [];
    C(i).Xend = [];
    C(i).X = [];
end

% Calculate average baseline activity in 10-D space
yB = [C.yB];
yB = mean(yB,2);

% Calculate average baseline activity in high-D space
zB = [C.zB];
zB = mean(zB,2);

% Loop over trials, calculate the distance from baseline
for i = 1:nTrials
    % Calculate distance from baseline for all samples
    y = C(i).y - repmat(yB,1,C(i).nSamp);
    z = C(i).z - repmat(zB,1,C(i).nSamp);
    d_y = nan(1,C(i).nSamp);
    d_z = nan(1,C(i).nSamp);
    for j = 1:C(i).nSamp
        d_y(j) = norm(y(:,j));
        d_z(j) = norm(z(:,j));
    end
    C(i).d_z = d_z;
    C(i).d_y = d_y;
    
    % Find time index with max distance from baseline.  Make sure this
    % search is done after cue/trajectory onset.
    d_z = d_z((C(i).onset(1) + 1):end);
    d_y = d_y((C(i).onset(1) + 1):end);
    [~,max_ind_z] = max(d_z);
    [~,max_ind_y] = max(d_y);
    C(i).max_ind_z = max_ind_z;
    C(i).max_ind_y = max_ind_y;
end

% Find time bin with most points.  This will be defined with respect to the
% appearance of the target.
edges = .5:100;
[N] = histcounts([C.max_ind_y],edges);
[~,max_ind_y] = max(N);

% Determine which data points to use for the movement data.  Loop over all
% trials and identifty the data depending on the 'kinOption' variable.
% This will add the kinematic data to the 'C' structure.
for i = 1:nTrials
    C(i) = getCalibrationData(C(i), kinOption, 'max_ind', max_ind_y,...
        'nStartSamp',nStartSamp,'nTargSamp',nTargSamp);
end

X = [C.X] * kinScaleFactor;
Y = [C.Y];
X = X(1:2,:); % Remove z-dimension

% Plot results

% Setup figure
n_row = 1;
n_col = 2;
ax_sz = 300;
ax_sp = 75;
[fw,fh,Ax] = calcFigureSize(n_row, n_col, ax_sz, ax_sz, ax_sp);
F(3) = figure('Position',[10 10 fw fh]);
set(F(3),'Name',[P.decoderName '_ModulationWindow'])

% Subplot 1: distance from baseline vs time
subplotSimple(n_row, n_col, 1, 'Ax', Ax); hold on;
for i = 1:nTrials
    plotDistanceFromBaseline(C(i))
end
plot([max_ind_y max_ind_y],get(gca,'YLim'),'--', ...
    'color',[0.75, 0, 0.75], ...
    'LineWidth', 2)
set(gca,'TickDir','out','box','off')
xlabel(sprintf('Bin (%d ms, aligned to cue onset)',P.binwidth), ...
    'FontSize', 12)
ylabel('Distance from baseline - latent space (a.u.)', ...
    'FontSize', 12)
title('Distance from baseline', 'FontSize', 14)

% Subplot 2: histogram of indices where distance from baseline was highest
subplotSimple(n_row, n_col, 2, 'Ax', Ax); hold on;
%H = histogram('BinEdges', edges, 'BinCounts', N);
%set(H,'FaceColor',ones(1,3)*.75)
plot([max_ind_y max_ind_y],get(gca,'YLim'),'--', ...
    'color',[0.75, 0, 0.75], ...
    'LineWidth', 2)
set(gca,'XLim',[0 round(1500/P.binwidth)],'TickDir','out','box','off')
xlabel(sprintf('Bin (%d ms, relative to cue onset)',P.binwidth), ...
     'FontSize', 12)
ylabel('Counts', 'FontSize', 12)
title('Max distance from baseline', 'FontSize', 14)

titleStr = sprintf('%s %s :: Decoder %d :: Calibration data modulation window', ...
    P.subject, P.date, P.decoderID);
plotTitle(titleStr)

fprintf('done.\n')

%--------------------------------------------------------------------------
% Remove fast latents (if desired)

tauThresh = G.GPFAParams.binWidth * tauBinThresh; % Threshold (ms)

% Convert gamma to time constant
tau = G.GPFAParams.binWidth ./ sqrt(G.GPFAParams.gamma); % time constant (ms)

% Remove latent dimensions with fast time constants -- Currently not used.
% Instead, the 'excludeLatents' variable can be used to remove specific
% dimensions if desired.
tauMask = tau >= tauThresh;

% Remove selected latents
xDim = size(Y,1);
tauMask = ~ismember(1:xDim,excludeLatents);
Y = Y(tauMask,:);

%--------------------------------------------------------------------------
% Calibrate decoder (OLE)

fprintf('Calibrating decoder ... ')

% Fit observation model based on data
% Note: might want to fit this better (make sure even number of each target)
yMean = mean(Y,2); 

xDim = size(Y,1);
kinDim = size(X,1);
N = size(Y,2);
Xfit = [ones(1,N) ; X];
B = nan(xDim,kinDim + 1);
pVals = nan(xDim,1);
for i = 1:xDim
    % Fit observation model for each dimension independently
    [b,~,~,~,STATS] = regress(Y(i,:)',Xfit');
    B(i,:) = b;
    pVals(i) = STATS(3);
end

% Find decoding weights as pseudoinverse of observation model
B = B(:,2:end);
W_ole = (inv(B'*B))*B';
c_ole = zeros(kinDim,1);

fprintf('done.\n')

%--------------------------------------------------------------------------
% Calibrate decoder (LR)

% Find baseline offset.  This is calculated as the average neural activity
Ystart = [C.Ystart];
Ystart = Ystart(tauMask,:);

d = mean(Ystart,2);
Y = Y - repmat(d,1,size(Y,2));

W_lr = X*(Y')*inv(Y*Y');
c_lr = -W_lr * d;

%--------------------------------------------------------------------------
% Select decoding weights to use

% Plot decoding weights
F(4) = figure();
set(F(4),'Name','DecoderWeightsComparison')

subplot(2,1,1)
imagesc(W_ole)
set(gca,'XTick',1:xDim,'YTick',1:2)
colormapNew(gca,'blue-orange');
set(gca,'CLim',[-25 25])
title('OLE')

subplot(2,1,2)
imagesc(W_lr)
set(gca,'XTick',1:xDim,'YTick',1:2)
colormapNew(gca,'blue-orange');
set(gca,'CLim',[-25 25])
xlabel('Latent state')
ylabel('Kinematic state')
title('Linear regression')

weightOption = 'lr';
switch weightOption
    case 'ole'
        W = W_ole;
        c = c_ole;
    case 'lr'
        W = W_lr;
        c = c_lr;
end

% Apply orthogonalization if necessary.  NOTE: if latents are to be
% excluded, the 'C' matrix here will probably have to be truncated
% appropriately.
switch xSpec
    case 'xorth'
        % Orthogonalize loading matrix to get transformation matrix from
        % non-orthonormalized to orthonormalized spaces
        C = G.GPFAParams.C;
        C = C(:,tauMask);
        [~,~,TT] = orthogonalize(zeros(G.GPFAParams.xDim,1),C);
        W = W * TT;
    case 'xsm'
        % Do nothing 
end

% Zero-pad W to full latent dimensionality
xDim = G.GPFAParams.xDim;
Wpad = zeros(kinDim,xDim);
Wpad(:,tauMask) = W;
W = Wpad;

P.W = W;
P.c = c;

%--------------------------------------------------------------------------
% Predict training data

fprintf('Predicting trajectores for calibration data ... ')

% Get RT GPFA parameters
RTParams = computeCausalRTParams(G.GPFAParams,'Tmax',P.Tmax,'nBins',P.nBins);

% Get GPFA parameters
P.nBins = RTParams.nBins;
P.xDim = RTParams.xDim;
P.M = RTParams.M;
P.CRinv = RTParams.CRinv;
P.d = RTParams.d;

% Predict cursor trajectories for calibration data
TDpredict = predictDecodeState_GPFA(G.TD,P);

% Find average trajectories (used for plotting)
TDpredictAvg = TDpredict.average;

fprintf('done.\n')

%--------------------------------------------------------------------------
% Put parameters into structure

% Zero-pad parameters up to full dimensionality.  This should be 96*2,
% where the order is [Ch1Srt1 Ch1Srt2 Ch2Srt1 Ch2Srt2 ... ].
nCh = 96;       % Hard-coded for now
nSort = 2;      % Determined by TDT capability
xDim = G.GPFAParams.xDim;
nUnits = nCh*nSort;
allSort = repmat((1:nSort)',nCh,1);
allCh = reshape(repmat((1:nCh)',1,nSort)',[],1);
chMask = ~ismember(allCh,P.exclCh);
sortMask = ismember(allSort,P.validSort);
unitMask = chMask & sortMask;

% Zero-pad parameters
CRinv = zeros(xDim,nUnits);
CRinv(:,unitMask) = RTParams.CRinv;
d = zeros(nUnits,1);
d(unitMask) = RTParams.d;

% Invert Y-axis.  Calibration data is inverted when loaded, so this needs
% to be un-done when going back to MonkeyHost.
W(2,:) = -W(2,:);
c(2) = -c(2);

P.M = RTParams.M;
P.CRinv = CRinv;
P.d = d;
P.W = W;
P.c = c + P.workspaceCenter';
P.G = G;

%--------------------------------------------------------------------------
% Plot results

fprintf('Plotting results ... ')

% Create figure
F(5) = figure('Position',[100 100 900 450]);

% Plot averaged GPFA trajectories
subplot(1,2,1)
TDpredictAvg.plot
title('GPFA Decoder Trajectories (averaged)')

% Plot all GPFA trajectories
subplot(1,2,2)
TDpredict.plot
title('GPFA Decoder Trajectories (all trials)')

% Add dataset info and set figure name
plotDatasetInfo('Earl',G.TD(1).date,[G.TD.trialID])
saveStr = sprintf('Earl_%s_GPFADecoderTrajectories', ...
    datestr(G.TD(1).date,'yyyymmdd'));
set(F(5),'Name',saveStr)

fprintf('done.\n')

end


function C = getCalibrationData(C, kinOption, varargin)
%
%

% Optional arguments
max_ind = [];  % Max. distance index ('maxDist' only)
nStartSamp = 3;  % Number of start samples to use
nTargSamp = 5;  % Number of target/movement samples to use

assignopts(who, varargin);

% Define kinematic data
switch kinOption
    case 'maxDist'
        % Maximum distance from baseline.  This case takes data surrounding
        % the time step when neural activity was most frequently furthest
        % from the baseline distribution.
        
        % Determine onset and offset indices.  Need to add the cue onset to
        % the 'max_ind' here, as 'max_ind' is defined with respect to cue
        % onset.
        onsetInd = max_ind + C.onset(1) - floor(nTargSamp/2);
        offsetInd = max_ind + C.onset(1) + floor(nTargSamp/2);
        if offsetInd > C.nSamp
            offsetInd = C.nSamp;
            onsetInd = offsetInd - nTargSamp + 1;
        end

    case 'endOfTrial'
        % End of trial.  This case takes samples from the end of the trial
        % -- right before the target was acquired.  In this case the offset
        % of the trial needs might need to be accounted for -- for Earl's
        % experiments this wasn't done, so this is being left alone for the
        % sake of consistency.
        
        startOnset = 1;
        startOffset = startOnset + nStartSamp - 1;
        
        offsetInd = C.nSamp;
        onsetInd = C.nSamp - nTargSamp + 1;

    case 'targetRamp'
        % Target ramp.  This was intended to allow all data to be used for
        % training.  THIS HAS BEEN DEPRECIATED AND IS NO LONGER FUNCTIONAL.
        
        % Throw error -- this case is no longer supported
        error('The "targetRamp" kinematic option is no longer supported.')
        
        centerPos = C.startPos;
        targPos = C.endPos;
        nCenterTs = 3;
        nRampTs = 3;

        % Generate linear ramp from center to target
        dPos = targPos - centerPos;
        dPosStep = dPos/nRampTs;
        dPosStep = repmat(dPosStep,1,nRampTs);
        rampPos = cumsum(dPosStep,2);

        pos = [repmat(centerPos,1,nCenterTs), ...
            rampPos, repmat(targPos,1,(C.nSamp - nCenterTs - nRampTs))];

        % Add data to structure
        C.X = pos;
        C.Y = C.y;
        C.Ystart = yStart;
        C.Xstart = repmat(C.startPos,1,nCenterTs);

    case 'initialPush'
        % Get data from the initial time steps after those used for the
        % center target.
        %startOnset = 1;
        %startOffset = nStartSamp;

        % Take neural activity centered about cue onset for baseline
        %startOnset = C.onset(1) - 1;
        %startOffset = startOnset + nStartSamp - 1;
        startOnset = C.onset(1);
        startOffset = startOnset + nStartSamp - 1;
        
        % Set onset and offset samples.  Make sure that the offset does
        % not exceed the number of total samples
        onsetInd = C.onset(1) + 3;  % Want to start 2 samples after cue onset
        offsetInd = onsetInd + nTargSamp - 1;
%         onsetInd = C.onset(1) + 8;
%         offsetInd = onsetInd + 2;
        
        if offsetInd > C.nSamp
            % Desired offset index exceeds the number of samples, truncate
            % to the total number of available samples.
            offsetInd = C.nSamp;
            
            % Check to see if the onset index is greater than the offset
            % index.  This can happen if the trial is successful really
            % quickly.  In this case, set the onset index to the onset of
            % the cursor control state.
            if onsetInd > offsetInd
                onsetInd = min(C.onset(2), offsetInd);
                warning('Desired offset for trial %d is before the specified onset.  The onset has been adjusted.', ...
                    C.trialID)
            end
        end
        
    case 'hybrid'
        % This approach is a hybrid of the 'initialPush', and 'endOfTarget'
        % cases.  It takes data from the end of the trial for successful
        % trials, and uses the initial push for failed trials.
        
        % Get data from the start of the trial up until ~150ms after cue
        % onset.  This is used to define the resting position/center of the
        % decoder.
        %startOnset = C.onset(1);
        %startOffset = C.onset(2) - 1;
        startOnset = C.onset(1);
        startOffset = startOnset + nStartSamp - 1;
        
        % Get target data.  Use different modes depending on whether or not
        % the trial was successful
        if C.successful
            % Trial was successful.  Get neural data from the end of the
            % trial.
            offsetInd = length(C.d_y);
            onsetInd = offsetInd - 3 + 1;
            
            % If desired onset is less than the onset of cursor control,
            % set it equal to the onset of the cursor control state
            onsetInd = max(C.onset(2), onsetInd);
        else
            % Trial was not successful.  Get data from the initial push
            % window.
        
            % Set onset and offset samples.  Make sure that the offset does
            % not exceed the number of total samples
            onsetInd = C.onset(1) + 8;
            offsetInd = onsetInd + nTargSamp - 1;

            if offsetInd > C.nSamp
                % Desired offset index exceeds the number of samples, truncate
                % to the total number of available samples.
                offsetInd = C.nSamp;

                % Check to see if the onset index is greater than the offset
                % index.  This can happen if the trial is successful really
                % quickly.  In this case, set the onset index to the onset of
                % the cursor control state.
                if onsetInd > offsetInd
                    onsetInd = min(C.onset(2), offsetInd);
                    warning('Desired offset for trial %d is before the specified onset.  The onset has been adjusted.', ...
                        C.trialID)
                end
            end
        end
        
end

% Get neural and kinematic data
yStart = C.y(:,startOnset:startOffset);
yEnd = C.y(:,onsetInd:offsetInd);

C.Xstart = repmat(C.startPos,1,size(yStart,2));
C.Xend = repmat(C.endPos,1,size(yEnd,2));
C.X = [C.Xstart C.Xend];

C.Ystart = yStart;
C.Yend = yEnd;
C.Y = [yStart yEnd];

% Define neural data mask for baseline
mask_base = false(1,C.nSamp);
mask_base(startOnset:startOffset) = true;
C.mask_base = mask_base;

% Define neural data mask for movement
mask_move = false(1,C.nSamp);
mask_move(onsetInd:offsetInd) = true;
C.mask_move = mask_move;
    
end


function plotDistanceFromBaseline(C)

% Get data to plot
d = C.d_y;

% Align to cue onset -- note that the 'max_ind' variable is aligned to
% cue onset, so need to add the onset back to this to get the approprite
% index into the distance variable.
onset = C.onset;
onset(onset > C.nSamp) = C.nSamp;  % Due to rounding this can get off slightly
x = (1:length(d)) - onset(1);
x_max_ind = C.max_ind_y;  % max_ind is already aligned
y_max_ind = d(x_max_ind + onset(1));
x_cue_onset = 0;
y_cue_onset = d(onset(1));
x_reach_onset = onset(2) - onset(1);
y_reach_onset = d(onset(2));

% Define masks
mask_base = C.mask_base;
mask_move = C.mask_move;

% Plot all data data (in gray)
plot(x,d,'color',ones(1,3)*.5)

% Plot baseline data (in blue)
plot(x(mask_base),d(mask_base),'color',[0, 0, 0.75], 'LineWidth', 2)

% Plot movement data (in purple)
plot(x(mask_move),d(mask_move),'color',[0.75, 0, 0.75], 'LineWidth', 2)

% Plot maximum distance from baseline (in purple)
plot(x_max_ind,y_max_ind, ...
    'color',[0.5, 0, 0.5], ...
    'Marker','.', ...
    'MarkerSize',15)
% Plot cue onset (black)
plot(x_cue_onset, y_cue_onset, ...
    'color', 'k', ...
    'Marker', '.', ...
    'MarkerSize', 10)
% Plot start of reach state (green)
plot(x_reach_onset, y_reach_onset, ...
    'color', [0, 0.75, 0], ...
    'Marker', '.', ...
    'MarkerSize', 10)

end
