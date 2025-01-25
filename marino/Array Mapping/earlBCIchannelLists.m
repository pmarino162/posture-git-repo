clear; clc; clf; close all;

%% Get maps for both flex and omnetics configs
    config = 'earlRH_LGAOmnetics_FlexConnector';
    [tdt2adapterMapFlex, adapter2cereportMapFlex, cereport2arrayMapFlex,array2anotomyMapFlex,tdt2arrayMapFlex,tdt2anatomyMapFlex] = getMaps(config);
    
%% Input anatomical sections in array pins
    PMd_array = 1:64;
    M1_array = 65:128;
    postPMd_array = [2:8,11:16,20:24,29:32,38:40,47:48,56];
    midPMd_array = [1 10 19 28 37 46 55 64];
    antPMd_array = [9,17:18,25:27,33:36,41:45,49:54,57:63];
    postAndMidPMd_array = union(postPMd_array,midPMd_array);
    postAndMidPMdInAnatomicalOrder_array = [8 7 16 6 15 24 5 14 23 32 4 13 22 31 40 3 12 21 30 39 48 2 11 20 29 38 47 56 1 10 19 28 37 46 55 64];
    
% %% Input TDT channels tuned to delay and reach
%     delay_tdt = [2:13,97:120];
%     reach_tdt = [65:78];
%     
%% Get anatomical sections in TDT channels
    %M1
    M1_tdt = nan(1,size(M1_array,2));
    for i = 1:length(M1_tdt)
       M1_tdt(i) = tdt2arrayMapFlex(tdt2arrayMapFlex(:,2)==M1_array(i),1); 
    end
    M1_tdt = sort(M1_tdt);
    %PMd
    PMd_tdt = nan(1,size(PMd_array,2));
    for i = 1:length(PMd_tdt)
       PMd_tdt(i) = tdt2arrayMapFlex(tdt2arrayMapFlex(:,2)==PMd_array(i),1); 
    end
    PMd_tdt = sort(PMd_tdt);
    %Post and Mid PMd in anatomical order (posterior->anterior)
    postAndMidPMdInAnatomicalOrder_tdt = nan(1,size(postAndMidPMdInAnatomicalOrder_array,2));
    for i = 1:length(postAndMidPMdInAnatomicalOrder_tdt)
       postAndMidPMdInAnatomicalOrder_tdt(i) = tdt2arrayMapFlex(tdt2arrayMapFlex(:,2)==postAndMidPMdInAnatomicalOrder_array(i),1); 
    end
    postAndMidPMdInAnatomicalOrder_tdt = sort(postAndMidPMdInAnatomicalOrder_tdt);
    %Post PMd
    postPMd_tdt = nan(1,size(postPMd_array,2));
    for i = 1:length(postPMd_tdt)
       postPMd_tdt(i) = tdt2arrayMapFlex(tdt2arrayMapFlex(:,2)==postPMd_array(i),1); 
    end
    postPMd_tdt = sort(postPMd_tdt);
    %Post and Mid PMd
    postAndMidPMd_tdt = nan(1,size(postAndMidPMd_array,2));
    for i = 1:length(postAndMidPMd_tdt)
       postAndMidPMd_tdt(i) = tdt2arrayMapFlex(tdt2arrayMapFlex(:,2)==postAndMidPMd_array(i),1); 
    end
    postAndMidPMd_tdt = sort(postAndMidPMd_tdt);
    
% %% Get delay and reach channels in array pins
%     delay_array = nan(1,size(delay_tdt,2));
%     reach_array = nan(1,size(reach_tdt,2));