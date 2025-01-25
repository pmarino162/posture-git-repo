function [saveDir] = getExperimentSaveDir20220419(dataset)

%This function takes dataset names in the format [monkey, date] (e.g.
%'E20200317')and outputs a string for the filepath for the folder in which
%the dataset is saved. 

monkey = dataset(1);
switch monkey
    case 'E'
        monkey = 'Earl';
    case 'N'
        monkey = 'Nigel';
    case 'R'
        monkey = 'Rocky';
end

year = dataset(2:5);

month = dataset(6:7);

date = dataset(2:9);

saveDir = fullfile('D:\Animals',monkey,year,month,date);

end