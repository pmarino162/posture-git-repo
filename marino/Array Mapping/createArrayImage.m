function [] = createArrayImage(config,chData,dateStr)

%% Load map
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap,array2anotomyMap,tdt2anatomyMap] = getMaps(config);

%% Noise channels
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
    
%% Get Clims and labels
    lowerLim = prctile(chData(chData>0),5);%min(chData(chData > 0));
    upperLim = prctile(chData(chData>0),95);%max(chData);
    clims = [lowerLim,upperLim];
    lowerLabel = round(ceil(clims(1,1)));
    upperLabel = round(floor(clims(1,2)));
    
%% Plotting 
    f = figure
    f.Position = [0 0 500 1000]
    cmap = parula;
    axis equal

    %PMd Plot
    subplot(2,1,1)
    imagesc(PMdArrayImage,clims)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
%     for row = 8:16
%         for col = 1:8
%             if ~isempty(find(
%                 text(j-0.1,i,num2str(anatomical2tdtMap(i,j)),'Fontsize',10,'Color',[1 1 1])
%             end
%         end
%     end
     c = colorbar('Ticks',[lowerLabel,upperLabel],'TickLabels',{['<=',num2str(lowerLabel)],['>=',num2str(upperLabel)]});
    title('PMd Array')
    
    %M1 Plot
    subplot(2,1,2)
    imagesc(M1ArrayImage,clims)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
%     for i=row:8
%         for col=1:8
%             if anatomical2tdtMap(i+8,j) ~= 0
%                 text(j-0.1,i,num2str(anatomical2tdtMap(i+8,j)),'Fontsize',10,'Color',[1 1 1])
%             end
%         end
%     end
      
      c = colorbar('Ticks',[lowerLabel,upperLabel],'TickLabels',{['<=',num2str(lowerLabel)],['>=',num2str(upperLabel)]});
%     c.Label.String = 'Vpp (uV)';
    title('M1 Array')
    sgtitle(dateStr)
end