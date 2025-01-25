clear;clc;clf;close all;

%% Load mappings 
% config = 'PMd';
% dateStr = '20210530';

config = 'M1';
dateStr = '20210530 M1';
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap, array2anatomyMap] = getMaps(config);
    
%% Get TDT Channel Quality
    tdtChQuality = ones(1,96);
    switch config
        case 'PMd'
        nonNeural = [23:25,30,31,38,55,59,64,67,71,77];    
        case 'M1'
        nonNeural = [9,14,15,24,27,33,48,61,67,79,80,85,93];
    end
    tdtChQuality(nonNeural) = 0;


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
    cmap = [1 0 0; .2 .5 .2];
    axis equal
    clims = [0,1];
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
    c.Label.String = 'Vpp (uV)';
    colormap(cmap)
    title('M1 Array')
    sgtitle(dateStr)
%     saveas(f,[savePath,'ArrayMapNeural',dateStr,'.jpg'])