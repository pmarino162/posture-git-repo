function [] = int2formatted_Marino(dataset)

%Loads intermediate Chase Lab Datasets and uses Chase Lab fxns to convert
%them to formatted data.  Then saves them. 


switch dataset
    case 'R20201020'
        %Load and convert data
        intData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390_intermediate.mat');
        intData = intData.Data;
        [Data, ~] = Intermediate2Formatted_v3(intData,[]);
        %Clean up rawSpikes
        [Data] = cleanSpikesRockyBC(Data);     
        %Save formatted data    
        save('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00390.mat','Data');
        
        %Load and convert data
        intData = load('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00393_intermediate.mat');
        intData = intData.Data;
        [Data, ~] = Intermediate2Formatted_v3(intData,[]);
        %Clean up rawSpikes
        [Data] = cleanSpikesRockyBC(Data);     
        %Save formatted data    
        save('D:\Animals\Rocky\2020\10\20201020\Rocky.VR.00393.mat','Data');
             
    case 'R20201021'
        %Load and convert data
        intData = load('D:\Animals\Rocky\2020\10\20201021\Rocky.VR.00394_intermediate.mat');
        intData = intData.Data;
        [Data, ~] = Intermediate2Formatted_v3(intData,[]);
        %Clean up rawSpikes
        [Data] = cleanSpikesRockyBC(Data);     
        %Save formatted data    
        save('D:\Animals\Rocky\2020\10\20201021\Rocky.VR.00394.mat','Data');
        
        %Load and convert data
        intData = load('D:\Animals\Rocky\2020\10\20201021\Rocky.VR.00396_intermediate.mat');
        intData = intData.Data;
        [Data, ~] = Intermediate2Formatted_v3(intData,[]);
        %Clean up rawSpikes
        [Data] = cleanSpikesRockyBC(Data);     
        %Save formatted data    
        save('D:\Animals\Rocky\2020\10\20201021\Rocky.VR.00396.mat','Data');
end




end