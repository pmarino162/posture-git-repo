clear; clc; close all;

   isoDatasetList = {'E20200116','E20200117','E20200120'};
   
   multipleTasksDatasetList = {'E20200314'};
   
   bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','E20210901','N20171215','N20180221','R20201020','R20201021'};
   
   reachDatasetList = {'E20210706','E20210707','E20210708','E20210709','E20210710',...
       'N20190222','N20190226','N20190227','N20190228','N20190301','N20190307',...
       'R20200221','R20200222'};
  
   controlDatasetList = {'E20190830','E20190729'};
   
   for datasetList = isoDatasetList
         dataset = datasetList{1,1}
         [Data,zScoreParams] = preprocessAndSaveData20220419(dataset);

   end