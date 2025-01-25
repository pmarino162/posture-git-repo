% clear;
% clc;
% clf;
% close all;

%% Load mappings 
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap, array2anatomyMap] = getMaps(config);
    
%% Get TDT Channel Quality
    tdtChQuality = zeros(1,96);
%     %Input data
%     goodChannels = [1:96];
% %      offChannels = [12,52];
%     offChannels = [];
%     lowAmpChannels = [5,7,10,18,19,22,25,26,27,32,36,37,42,44,51,54,60,61,62,65,...
%         71,73,75,77,79,90,91,92,94];   
%     intermediateCh = [];
%     badCh = [5,10,12,14:16,18,20:24,];
%   add color for not tested  
%      tdtChQuality = loadChQuality('20210524_channelQuality.txt');
%     tdtChQuality = loadChQuality('20210525_channelQuality.txt');
%     tdtChQuality25 = loadChQuality('20210525_channelQuality.txt');
tdtChQuality = Vpp;
switch dateStr
    case '20210524 M1'
        tdtChQuality([5,11:12,61:76,81:83,93:96]) = -10;
    case '20210528 M1'
        tdtChQuality([10 12 20 30 52]) = -10;
    otherwise
        tdtChQuality([12 52]) = -10;
end



%% Map from TDT to electrodes    
    arrayChQuality = zeros(1,128);
    tdt2arrayMap = zeros(96,2);
    for i = 1:96
        adapterInd = find(cell2mat(tdt2adapterMap(:,1))==i);
        adapterPin = tdt2adapterMap{adapterInd,2};
        cereportInd = find(contains(adapter2cereportMap(:,1),adapterPin));
        cereportPin = adapter2cereportMap{cereportInd,2};
        arrayInd = find(contains(cereport2arrayMap(:,1),cereportPin));
        arrayPin = cereport2arrayMap{arrayInd,2};
        tdt2arrayMap(i,1) = i; tdt2arrayMap(i,2) = arrayPin;
        
        arrayChQuality(arrayPin) = tdtChQuality(i);

    end



%% Create array images
anatomicalImage = zeros(16,8);

for row = 1:16
    for column = 1:8  
        arrayPin = array2anatomyMap(row,column);
        anatomicalImage(row,column) = arrayChQuality(arrayPin);
    end
end

PMdArrayImage = anatomicalImage(1:8,:);
M1ArrayImage = anatomicalImage(9:16,:);

%% Define anatomical2tdtMap
anatomical2tdtMap = zeros(16,8);
for row = 1:16
    for column = 1:8 
        arrayPin = array2anatomyMap(row,column);
        if ismember(arrayPin,tdt2arrayMap(:,2))
            tdtCh = find(tdt2arrayMap(:,2)==arrayPin);
            anatomical2tdtMap(row,column) = tdtCh;
        end
    end
end

%% Plotting 
savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Choking\Prez Array Characterization';
    f = figure
    f.Position = [0 0 500 1000]
%     cmap = zeros(65,3);
%     cmap(2:65,:) = parula;
%     cmap(1,:) = [1 0 0];
    cmap = parula;
    axis equal
    clims = [40,400];
    %PMd Plot
    subplot(2,1,1)
    imagesc(PMdArrayImage,clims)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
    yticks = [];
    for i=1:8
        for j=1:8
            if anatomical2tdtMap(i,j) ~=0
            text(j-0.1,i,num2str(anatomical2tdtMap(i,j)),'Fontsize',10,'Color',[1 1 1])
            end
        end
    end
%     c = colorbar('Ticks',[-10,10,round(clims(1,2)/2),clims(1,2)],'TickLabels',{'Noise','Off (0)',num2str(round(clims(1,2)/2)),num2str(round(clims(1,2)))});
    c = colorbar('Ticks',[clims(1,1),clims(1,2)],'TickLabels',{'<=40','>=400'});
    c.Label.String = 'Vpp (uV)';
    colormap(cmap)
    title('PMd Array')
    %M1 Plot
    subplot(2,1,2)
    imagesc(M1ArrayImage,clims)
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
    yticks = [];
    for i=1:8
        for j=1:8
            if anatomical2tdtMap(i+8,j) ~= 0
            text(j-0.1,i,num2str(anatomical2tdtMap(i+8,j)),'Fontsize',10,'Color',[1 1 1])
            end
        end
    end
%     c = colorbar('Ticks',[-10,10,round(clims(1,2)/2),clims(1,2)],'TickLabels',{'Noise','Off (0)',num2str(round(clims(1,2)/2)),num2str(round(clims(1,2)))});
    c = colorbar('Ticks',[clims(1,1),clims(1,2)],'TickLabels',{'<=40','>=400'});
    c.Label.String = 'Vpp (uV)';
    colormap(cmap)
    title('M1 Array')
    sgtitle(dateStr)
    saveas(f,[savePath,'ArrayMap',dateStr,'.jpg'])