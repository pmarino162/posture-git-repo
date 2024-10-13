% Define the directories to search
dirs = {
    'D:\Animals\Earl\2020\01\20200116', ...
    'D:\Animals\Earl\2020\01\20200117', ...
    'D:\Animals\Earl\2020\01\20200120'
};

% Define the folder name patterns to search for
validPatterns = {'I15', 'I30', 'N', 'E15', 'E30'};

% Initialize an array to store all scaling values
allScalingValues = [];
allTargetYValues = [];
% Loop through each directory
for i = 1:length(dirs)
    currentDir = dirs{i};
    
    % Get a list of all subfolders
    folderInfo = dir(currentDir);
    
    % Loop through each folder
    for j = 1:length(folderInfo)
        folderName = folderInfo(j).name

        % Check if the folder name matches any of the valid patterns
        if folderInfo(j).isdir && any(endsWith(folderName, validPatterns))
            folderPath = fullfile(currentDir, folderName);
            
            % Find the 'translated.mat' file in the folder
            matFile = dir(fullfile(folderPath, '*translated.mat'));
            
            if ~isempty(matFile)
                matFilePath = fullfile(folderPath, matFile(1).name);
                
                % Load the .mat file, which will create a variable named 'Data'
                load(matFilePath, 'Data');
                
                % Extract unique scaling values from Data.Definitions.forceTransformation.scaling(2)
                for row = 1:size(Data,2)
                    if isfield(Data(row).Definitions, 'forceTransformation') && ...
                       isfield(Data(row).Definitions.forceTransformation, 'scaling')
                       
                        % Get the scaling value
                        scalingValue = Data(row).Definitions.forceTransformation.scaling(2);
                        
                        % Append the scaling value to the list
                        allScalingValues = [allScalingValues, scalingValue];
                    end
                    
                    targetY = Data(row).Parameters.StateTable(3).StateTargets.location(1,2);
                    allTargetYValues = [allTargetYValues,targetY];
                end
            end
        end
    end
end

% Get unique scaling values
uniqueScalingValues = unique(allScalingValues);
uniqueTargetYValues = unique(allTargetYValues);
% Display the unique scaling values
disp('Unique scaling values:');
disp(uniqueScalingValues);
disp('Unique target values:')
disp(uniqueTargetYValues);