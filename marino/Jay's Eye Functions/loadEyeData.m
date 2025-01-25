function E = loadEyeData(D, dataDir)
% D can be struct or datestr (used to find filename)
% if D provided as struct, then adds trial and block indices for aligning

    if isstruct(D)
        dtstr = D.datestr;
    else
        dtstr = D;
        D = [];
    end
    if nargin < 2
        dataDir = '~/code/wmpCode/data/eyetracking/processed';
    end
    
    % find processed eye data files
    % may be multiple files (e.g. 'onlineDecode2' and 'onlineDecode3')
    fnms = dir(fullfile(dataDir, ['*' dtstr '*.mat']));
    if numel(fnms) == 0
        error('Could not find processed eye data.');
    end
    E = [];    
    for ii = 1:numel(fnms)
        infile = fullfile(dataDir, fnms(ii).name);
        e = load(infile);
        if isempty(E)
            E = e;
        else
            E.Data = [E.Data e.Data];
        end
    end
    if ~isempty(D)
        E.Data = findMatchingTrials(E.Data, D);
    end
    E = removeBlinks(E);
    E = addClockTime(E);
    E.datestr = D.datestr;
end

function E = addClockTime(E)
% given a datetime string, convert this to seconds, rel. to start of day
    tms = num2str(cellfun(@(c) c.date, {E.Data.Overview})');
    secs = tools.trialTimesToSeconds(tms);
    for ii = 1:numel(E.Data)
        E.Data(ii).trialTime = secs(ii);
    end
end

function E = removeBlinks(E, eyeName, pthresh, threshborder)
% 
% set analogData to nan whenever blink is detected
% 
    if nargin < 2
        eyeName = 'Left Pupil'; % this one tends to be more reliable
    end
    if nargin < 3
        pthresh = -9; % chosen by inspection to identify blink cluster
    end
    if nargin < 4
        threshborder = 40; % remove eye data 40 ms before and after blink
    end

    cind = ismember(E.dataNames, eyeName);
    for ii = 1:numel(E.Data)
        pup = E.Data(ii).TrialData.analogData(:,cind);
        n = numel(pup);
        blinks = find(pup < pthresh);
        
        for t = 1:numel(blinks)
            E.Data(ii).TrialData.analogData(...
                max(1,blinks(t)-threshborder):...
                min(n,blinks(t)+threshborder),:) = nan;
        end
    end    
end

function Data = findMatchingTrials(Data, D)
% Data and D may not necessarily contain the same trials
% this function finds the nearest trial_index from D for each trial in Data

    tmsData = nan(numel(Data),1);
    for ii = 1:numel(Data)
        tmsData(ii) = Data(ii).Overview.date;
    end
    tmsData = tools.trialTimesToSeconds(tmsData);
    
    tmsD = [];
    for jj = 1:numel(D.blocks)
        B = D.blocks(jj);
        blk = repmat(B.block_index, numel(B.trialTime), 1);
        trs = sort([unique(B.trial_index); B.trialsSkipped]);
        ctrs = unique(B.trial_index);
        ix = ismember(trs, ctrs);
        tmsD = [tmsD; B.trialTime(ix) ctrs blk(ix)];
    end
    ds = pdist2(tmsData, tmsD(:,1));
    [vs,ix] = min(ds,[],2);
    for ii = 1:numel(vs)
        tr = tmsD(ix(ii),:);
        if floor(vs(ii)) > 0
            % didn't find a close enough match;
            % this must be data not included in simpleData
            Data(ii).block_index = nan;
            Data(ii).trial_index = nan;
        else
            Data(ii).block_index = tr(3);
            Data(ii).trial_index = tr(2);
        end
    end
end

