
doSave = false;
dataDir = '~/code/wmpCode/data/eyetracking/raw';
saveDir = '~/code/wmpCode/data/eyetracking/processed';
fnms = dir(fullfile(dataDir, '*.mat'));

dataNamesToKeep = {'Left Eye X', 'Left Eye Y', ...
    'Right Eye X', 'Right Eye Y', 'Left Pupil', 'Right Pupil'};

if numel(fnms) > 1 && ~doSave
    error('Did you mean to not save?');
end
for ii = 1:numel(fnms)
    infile = fullfile(dataDir, fnms(ii).name);
    d = load(infile);
    if isfield(d, 'rawData')
        Data = d.rawData;
    else
        Data = d.Data;
    end
    dataNames = Data(1).Definitions.analogChannelNames;
    stateNames = cellfun(@(c) c{1}, ...
        {Data(1).Parameters.StateTable.stateName}, 'uni', 0);
    
    ix = ismember(dataNames, dataNamesToKeep);
    if sum(ix) ~= numel(dataNamesToKeep)
        warning(['Not all data fields found for ' fnms(ii).name]);
    end
    dataNames = dataNames(ix);
    
    datas = [];
    for jj = 1:numel(Data)
        data = Data(jj).TrialData;
        data = rmfield(data, 'Marker');
        data = rmfield(data, 'Decoder');
        data = rmfield(data, 'TDT');
        data.analogData = data.analogData(:,ix);
        Data(jj).TrialData = data;
    end
    Data = rmfield(Data, 'Definitions');
    Data = rmfield(Data, 'Parameters');    
    
    outfile = fullfile(saveDir, fnms(ii).name);
    save(outfile, 'Data', 'dataNames', 'stateNames');
end

