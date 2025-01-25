function saveGPFADecodingParams(P)

% Check to see if parameter saving would overwrite existing parameters
paramsFilename = fullfile(P.analysisDir,P.decoderName);
paramsFilename = [paramsFilename '.mat'];
dataFilename = ['MkHstCalibData' P.decoderName(9:end) '.mat'];
dataFilename = fullfile(P.analysisDir,dataFilename);

% Save parameters
if exist(dataFilename,'file') ~= 2
    fprintf('Saving parameters ... ')
    
    % Get data/parameters to save
    G = P.G;        % Get calibration data
    P.G = [];       % Remove calibration data from parameters
    bci_params = P; % Rename parameter structure for MonkeyHost
    
    % Save parameters/data
    save(paramsFilename,'bci_params')
    save(dataFilename,'G')
    
    fprintf('done.\n')
else
    warning('Save filename already exists.  Decoding parameters not saved.')
end