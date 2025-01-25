function [Data] = calibrateEyePosition(Data,varargin)
% CALIBRATEEYEPOSITION Takes eye voltage and translates it into eye position in Phasespace 
%           units using calibration matrix
%   [Data] = calibrateEyePosition(Data,horzgain,vertgain,horzoffset,vertoffset)
%
% 30-Jan-2014 BA: Commented code, renamed things to be more sensible, and
% optimized code. Has potential to be optimized more with proper use
% of parfor or some clever techniques. Compared with original calib_eye
% function. Data structure produced was exactly the same (verified with
% structeq in Chiou/usefulfunctions)

%% AUTHOR    : Erin Crowder, Berook Alemayehu, and others
%% $DATE     : 24-Jan-2014 $ 
%% $Revision : 2.00 $ 
%% DEVELOPED : 8.1.0.604 (R2013a) 

% ---------------------------BEGIN Input Parsing---------------------------
% The following 2 blocks of code enables the use of parameters, checks
% their values, and sets default values.
% You may change default values here.
p = inputParser;
p.addRequired('Data',@isstruct);
p.addParamValue('eye','Left',@(x) strcmpi(x,'left') | strcmpi(x,'right'));
p.parse(Data,varargin{:});

Data = p.Results.Data;
eyeSide = p.Results.eye;
% ---------------------------END Input Parsing---------------------------
removeTrials = [];

    for i = 1:length(Data)
        if Data(i).Overview.tag >= 0 % Setting the calibration matrix happens during tag 0.
            eyeVoltage = [];
            eyeSideIndex = find(strcmpi(Data(i).Definitions.analogChannelNames(:,1), [eyeSide ' Eye X']));
            eyeVoltageLength = length(Data(i).TrialData.analogData(:,1));
            eyeVoltage(1:eyeVoltageLength,2:3) = Data(i).TrialData.analogData(:,eyeSideIndex:eyeSideIndex+1); %altering analog eye matrix

            % For now, we use the fact that both eyes are moving the same to
            % deal w/ data in which only one eye was tracked by using the L eye
            % twice (when we have stereo data we can fix this & use both eyes)
            eyeVoltage(1:eyeVoltageLength,1) = 1;

            [aboveThresholdTimes,~] = find(abs(eyeVoltage(:,2:3)) > 9.3); %find blinks, voltage higher than 9.3 (voltage when eye is lost/not being tracked)
            aboveThresholdTimes = unique(aboveThresholdTimes);

            % If there are blinks in this trial, remove them as well as +/- 50ms
            if ~isempty(aboveThresholdTimes)
                blink_start = find(diff(aboveThresholdTimes)>1)';

                blinkTimeIndex = [1 blink_start + 1 length(aboveThresholdTimes)]; 
                % Since there may be more than 1 blink per trial, break up
                % blinks & handle separately
                for j=1:length(blinkTimeIndex)-1
                    
                    if blinkTimeIndex(j) == blinkTimeIndex(j+1)
                        ind = blinkTimeIndex(j);
                    else
                        ind = find(aboveThresholdTimes >= aboveThresholdTimes(blinkTimeIndex(j)) & aboveThresholdTimes < aboveThresholdTimes(blinkTimeIndex(j+1)));
                    end
                    
                    % Checks that there is 50ms before & after blink. Fills
                    % in with NaNs.
                    if aboveThresholdTimes(min(ind)) - 50 <= 0  & aboveThresholdTimes(max(ind)) + 50 > eyeVoltageLength
                        eyeVoltage((1:end),:) = NaN;
                    elseif aboveThresholdTimes(min(ind)) - 50 <= 0 
                        eyeVoltage((1 : aboveThresholdTimes(max(ind)) + 50),:) = NaN;
                    elseif aboveThresholdTimes(max(ind)) + 50 > eyeVoltageLength
                        eyeVoltage((aboveThresholdTimes(min(ind)) - 50 : end),:) = NaN;
                    else
                        eyeVoltage((aboveThresholdTimes(min(ind)) - 50 : aboveThresholdTimes(max(ind)) + 50),:) = NaN;
                    end
                end
            end

            % Copies monocular eye voltage and uses it as the second eye.
            % Calbration matrix assumes both eyes are moving the same.
            eyeVoltage(1:eyeVoltageLength,4:5) = eyeVoltage(1:eyeVoltageLength,2:3);
            
            calibMatrix = Data(i).Parameters.eyeTrackerCalibrationMatrix(:,:);
            if isempty(calibMatrix)
%                 Data(i).TrialData.EyePosition.raw = NaN;
%                 error('calibrateEyePosition:calibrationMatrix', ['No calibration matrix present during trial ' num2str(i)])
                fprintf(['No calibration matrix present during trial ' num2str(i) '. Trial ' num2str(i) ' will be removed.\n']);
                removeTrials = [removeTrials i];
            elseif sum(isnan(eyeVoltage)) == length(eyeVoltage)
                fprintf(['All eye voltages are above threshold during trial ' num2str(i) '. Removing trial ' num2str(i) ' will be removed.\n']);
                removeTrials = [removeTrials i];
            else
                raw = (calibMatrix*(eyeVoltage)')';
                Data(i).TrialData.EyePosition.raw = raw;
                
                if sum(isnan(raw(:,1))) > size(raw,1)*0.5
                    fprintf(['Warning: More than 50% of eye voltages were lost in trial: ', num2str(i), '\n']);
                end
            end



    %             end
    % 
    %             raw = Data(i).TrialData.EyePosition.raw;
    %             target_names = Data(i).Parameters.TrialTargets.names;
    %             windows = Data(i).Parameters.TrialTargets.window;
    %             start = find(ismember(target_names, 'start')==1,1,'first');
    %             reach = find(ismember(target_names, 'reach3')==1,1,'first');
    %             figure;
    %             plot(raw(:,1),-raw(:,2),'r.')
    %             rectangle('Position',[windows(start,1)-(0.5*windows(start,4)),-windows(start,2)-(0.5*windows(start,5)),...
    %                 windows(start,4),windows(start,5)])
    %             rectangle('Position',[windows(reach,1)-(0.5*windows(reach,4)),-windows(reach,2)-(0.5*windows(reach,4)),...
    %                 windows(reach,4),windows(reach,4)])
    %             daspect([1,1,1])
    %             
    %             eyepos = [eyepos; raw(:,1) -raw(:,2)];
    %             
    %             
    %             date = num2str(Data(1).Overview.date);
    %             filepath = fullfile('results',date(1:8),'behavior');
    %                 if ~isdir(filepath)
    %                     mkdir(filepath) 
    %                 end
    %             print('-dpdf',fullfile(filepath,['Eye Position Trial' num2str(i)]))
    %             close all
        end
    end
    
    Data(removeTrials) = [];
end