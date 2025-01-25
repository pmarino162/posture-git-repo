function [Data] = cleanSpikesRockyBC(Data)

%Cleans up raw spikes format imported from Chase Lab Intermediate Data
%format

        %Unpack spike info
        timestamps = Data.rawSpikes.timestamps; 
        channel = Data.rawSpikes.channel; 
        unit = Data.rawSpikes.unit; 
        source_index = Data.rawSpikes.source_index;
        
        %Remove timestamps for which time = 0 & channel = 0
        rmInd = timestamps==0 & channel ==0;
        timestamps(rmInd) = []; 
        channel(rmInd) = []; 
        unit(rmInd) = []; 
        source_index(rmInd) = []; 
        
        %Sort timestamps
        [timestamps,sortedInd] = sort(timestamps);
        channel = channel(sortedInd); 
        unit = unit(sortedInd); 
        source_index = source_index(sortedInd);
        
        %Align Timestamps
        
        %Save
        Data.rawSpikes.timestamps = timestamps; 
        Data.rawSpikes.channel = channel; 
        Data.rawSpikes.unit = unit; 
        Data.rawSpikes.source_index = source_index;

end