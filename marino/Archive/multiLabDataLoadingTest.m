clear; clc;
    
%% Compare preprocessed data from each script
    %Batista Earl Dataset
    dataset = 'E20200317';
%     getSpikes = true;
%         getSorts = false;
%         exclCh =  [3 12 16 20 31 73 85 92 94];
%         exclZero = true; %Exclude channels that only contain sort 0
%     getMarker = true;
%         checkPhasespaceSync = false;
%         droppedMarkerTrials = 'all';
%         removeBadPhasespace = false;
%         centerMarker = true; %Subracts workspace center from marker pose if true
%     getForce = true;
%         forceSetup = 'EL';
%     getAlg = false;
%     getEye = false;
%     getWaveforms = false;
%     getLFP = false;
%     getBB = false;
%     getKin = true;
%     checkDecode = true;
%     removeScreenFreezeTrials = true;
%     inclStateTable = false;
%     exclTrials = [];
% 
%     I30 = load('D:\Animals\Earl\2020\03\20200317\05_I30_grid\Earl20200317_05_I30_grid_SI_translated.mat');
%     rawData = I30.Data;
%     procData = preprocessBatistaData(rawData,'checkPhasespaceSync',checkPhasespaceSync,'droppedMarkerTrials',droppedMarkerTrials,'checkDecode',checkDecode,'removeScreenFreezeTrials',removeScreenFreezeTrials,'removeBadPhasespace',removeBadPhasespace);
    batistaData = loadData(dataset);
    %Chase Nigel BCI data
    [chaseData] = loadNigelData20180221_20220419;
    

%% Get number of trials for each of 2 datasets
%     datasets = {'N20180221','E20200317'};
    datasets = {'E20200317'};
    numDatasets = size(datasets,2);
    numTrials = zeros(1,numDatasets); 

    for setInd = 1:length(datasets)
         curSet = datasets{setInd};
         Data = loadData(curSet);
         numTrials(setInd) = size(Data,2);
    end