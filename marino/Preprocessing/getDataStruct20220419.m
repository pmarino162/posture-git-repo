function [Data] = getDataStruct20220419(rawData,dataset,varargin)

% This function processes raw experimental data so that it's in a
% user-friendly format

%% Variable Arguments
    getSpikes = false;
        getSorts = false;
        exclCh = [];
        exclZero = true; %Exclude channels that only contain sort 0
    getMarker = false;
        checkPhasespaceSync = false;
        droppedMarkerTrials = 'all';
        removeBadPhasespace = true;
        centerMarker = true; %Subracts workspace center from marker pose if true
    getForce = false;
        forceSetup = '';
        centerForceCursor = true; %Subracts workspace center from force cursor pose if true
    getAlg = false;
    getEye = false;
    getWaveforms = false;
    getLFP = false;
    getBB = false;
    getKin = false;
    checkDecode = false;
    removeScreenFreezeTrials = false;
    inclStateTable = false;
    exclTrials = [];
    assignopts(who,varargin);
    
%% Preprocess Data
    switch dataset
        %Batista Lab Data
        case {'E20200317'}
            procData = preprocessBatistaData(rawData,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'removeBadPhasespace',removeBadPhasespace);
            procData = removeCatchTrials(procData);
        %Chase Lab Data (Nigel)
        case {'N20180221'}
            procData = preprocessNigelBCIData(rawData);
    end
    %Clear rawData
    clearvars rawData
    
    
end