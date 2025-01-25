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
    for i = 1:length(Data)
        if Data(i).Overview.tag >= 0 % Setting the calibration matrix happens during tag 0.
            eyeVoltage = [];
            eyeSideIndex = find(strcmpi(Data(i).Definitions.analogChannelNames(:,1), [eyeSide ' Eye X']));
            eyeVoltageLength = length(Data(i).TrialData.analogData(:,1));
%             eyeVoltage(1:eyeVoltageLength,2:3) = Data(i).TrialData.analogData(:,[eyeSideIndex eyeSideIndex+1]); %altering analog eye matrix
            % NOT IN USE -> % HACK 3rd channel tends to be better Left Eye Y...
            eyeVoltage(1:eyeVoltageLength,2:3) = Data(i).TrialData.analogData(:,[eyeSideIndex eyeSideIndex+1]); %altering analog eye matrix

            % For now, we use the fact that both eyes are moving the same to
            % deal w/ data in which only one eye was tracked by using the L eye
            % twice (when we have stereo data we can fix this & use both eyes)
            eyeVoltage(1:eyeVoltageLength,1) = 1;

            [aboveThresholdTimes,~] = find(abs(eyeVoltage(:,2:3)) > 9.3); %find blinks, voltage higher than 9.3 (voltage when eye is lost/not being tracked)
            aboveThresholdTimes = unique(aboveThresholdTimes);

            % get the calibration matrix if present
            calibMatrix = Data(i).Parameters.eyeTrackerCalibrationMatrix(:,:);
            
            % If there are blinks in this trial, remove them as well as +/- 50ms
            if ~isempty(aboveThresholdTimes) & ~isempty(calibMatrix)
                % AS: It's only really a blink if you lose at least quite a
                % few in a row...otherwise they're random drops. Sigh.
                % We'll say you need 10ms in a row for it to be a blink
%                 blink_start = find(diff(aboveThresholdTimes)>1); % old way...
                blinkConsecIndCount = 10;
                blink_start = intersect(find([1; diff(aboveThresholdTimes)>1])',...         % same as before - when there hasn't been a drop in a while...
                    strfind([diff(aboveThresholdTimes)' 0],ones(blinkConsecIndCount,1)'));  % and there's at least 10 in a row. interesting that strfind works like this!
                if ~isempty(blink_start)
%                     blinkTimeIndex = [1 blink_start + 1 length(aboveThresholdTimes)];
                    blinkTimeIndex = blink_start;
                    % Since there may be more than 1 blink per trial, break up
                    % blinks & handle separately
                    % AS: Honestly I have no clue what was going on here,
                    % except for that any time there was a blink a trial
                    % was basically evicerated. Let's not do that.
                    for j=1:length(blinkTimeIndex)-1
                        
%                         if blinkTimeIndex(j) == blinkTimeIndex(j+1)
                            ind = blinkTimeIndex(j);
%                         else
%                             ind = find(aboveThresholdTimes >= aboveThresholdTimes(blinkTimeIndex(j)) & aboveThresholdTimes < aboveThresholdTimes(blinkTimeIndex(j+1)));
%                         end
                        
                        % Checks that there is 50ms before & after blink. Fills
                        % in with NaNs.
                        % AS: this is dumb. You're burning 50ms before a
                        % blink that's totally unnecessary. Let's not. Use
                        % 5ms instead.
                        if aboveThresholdTimes(min(ind)) - 5 <= 0  & aboveThresholdTimes(max(ind)) + 50 > eyeVoltageLength
                            eyeVoltage((1:end),:) = NaN;
                        elseif aboveThresholdTimes(min(ind)) - 5 <= 0
                            eyeVoltage((1 : aboveThresholdTimes(max(ind)) + 50),:) = NaN;
                        elseif aboveThresholdTimes(max(ind)) + 50 > eyeVoltageLength
                            eyeVoltage((aboveThresholdTimes(min(ind)) - 5 : end),:) = NaN;
                        else
                            eyeVoltage((aboveThresholdTimes(min(ind)) - 5 : aboveThresholdTimes(max(ind)) + 50),:) = NaN;
                        end
                    end
                end
                % AS: Any other above threshold, call nans.
                eyeVoltage(aboveThresholdTimes,:) = nan;
            end

            % Copies monocular eye voltage and uses it as the second eye.
            % Calbration matrix assumes both eyes are moving the same.
            % AS 20200519 no clue why they did this given the calibration
            % matrix is 3x3. W/e man...Also the calibration matrix is 3x3
            % implying a Z component..? Can you infer that from only two
            % sources...?
%             eyeVoltage(1:eyeVoltageLength,4:5) = eyeVoltage(1:eyeVoltageLength,2:3);
            
            if isempty(calibMatrix)
                Data(i).TrialData.EyeData.position = NaN;
%                 error('calibrateEyePosition:calibrationMatrix', ['No calibration matrix present during trial ' num2str(i)])
                fprintf(['No calibration matrix present during trial ' num2str(i) '. Trial ' num2str(i) ' will be nans.\n']);
                
            elseif sum(isnan(eyeVoltage(:,1))) == length(eyeVoltage)
                fprintf(['All eye voltages are above threshold during trial ' num2str(i) '. Removing trial ' num2str(i) ' will be removed.\n']);
                Data(i).TrialData.EyeData.position = NaN;
            else
                raw = (calibMatrix*(eyeVoltage)')';
                Data(i).TrialData.EyeData.position = raw;
                
                if sum(isnan(raw(:,1))) > size(raw,1)*0.5
                    fprintf(['Warning: More than 50 percent of eye voltages were lost in trial: ', num2str(i), '. Removing trial. \n']);
                    Data(i).TrialData.EyeData.position = nan;
                end
            end
        end
    end
end