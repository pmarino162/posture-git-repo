clear; clc; clf; close all;

%% Load data
    dateStr = '20210622';
    config = 'earlRHM1';
    chData = load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210624\Vpp20210622.mat');
    chData = chData.Vpp;
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap,array2anotomyMap,tdt2anatomyMap] = getMaps(config);

    chData(12) = 0;
    chData(52) = 0;
    
%% Create array images
    anatomicalImage = zeros(16,8);
    for channel = 1:size(tdt2anatomyMap,1)
       row =  tdt2anatomyMap{channel,2};
       col =  tdt2anatomyMap{channel,3};
       anatomicalImage(row,col) = chData(channel);  
    end
    PMdArrayImage = anatomicalImage(1:8,:);
    M1ArrayImage = anatomicalImage(9:16,:);
    rowList = cell2mat(tdt2anatomyMap(:,2));
    colList = cell2mat(tdt2anatomyMap(:,3));
    
%% Plotting 
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Choking\Earl Array Characterization';
    f = figure
    f.Position = [0 0 500 1000]
    cmap = parula;
    axis equal
    lowerLim = min(chData(chData > 0));
    upperLim = max(chData);
    clims = [lowerLim,upperLim];
    %PMd Plot
    subplot(2,1,1)
    imagesc(PMdArrayImage,clims)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
%     yticks = [];
%     for row = 8:16
%         for col = 1:8
%             if ~isempty(find(
%                 text(j-0.1,i,num2str(anatomical2tdtMap(i,j)),'Fontsize',10,'Color',[1 1 1])
%             end
%         end
%     end
     c = colorbar('Ticks',[round(clims(1,1)),round(clims(1,2))]);
%     c.Label.String = 'Vpp (uV)';
    title('PMd Array')
    
    %M1 Plot
    subplot(2,1,2)
    imagesc(M1ArrayImage,clims)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
%     yticks = [];
%     for i=row:8
%         for col=1:8
%             if anatomical2tdtMap(i+8,j) ~= 0
%                 text(j-0.1,i,num2str(anatomical2tdtMap(i+8,j)),'Fontsize',10,'Color',[1 1 1])
%             end
%         end
%     end
      c = colorbar('Ticks',[round(clims(1,1)),round(clims(1,2))]);
%     c.Label.String = 'Vpp (uV)';
    title('M1 Array')
    sgtitle(dateStr)
%     saveas(f,[savePath,'ArrayMap',dateStr,'.jpg'])