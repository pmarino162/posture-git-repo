function [Data] = findFixations(Data, varargin)
% FIND_FIX Finds fixation start and end points per trial.
%   Load preprocessed data first!
%   Data = find_fix(Data) takes preprocessed data, finds fixation start and
%       end points per trial, and updates the data
%   find_fix(...,'Parameter',ParameterValue,...) finds fixations with
%       available parameters listed below.
%
%   Parameters
%       showPlot - Boolean that if true, will show plots.
%               Default: 0 (false)
%       savePath - When set, saves figures to the specified path. Make sure
%               to write the '/' at the end.
%               Default: 'savePath',['Plots/' data.fileName '/']
%       saveFigs - Set true or false (or 1 or 0). 1: saves figures in
%           save path. 0: doesn't save figures.
%               Default: 1
%       mmConvert - Conversion factor from raw eye positions to
%           millimeters.
%               Default: 0.4
%       degConvert - Conversion factor from millimeter eye position to
%           degrees, in mm. Corresponds to distance from eye tracker to
%           subjects' eyes.
%               Default: 330

%   Examples:
%       Data = find_fix(Data,'showPlot',1,'saveFigs',1,'mmConvert',0.5)
%           finds fixation points, shows the plots and saves them in the
%           default save path, uses 0.5 conversion for mm instead of
%           default 0.4
%
%   NOTE: y-axis is flipped. Phasespace coordinate system reverses the
%       y-axis direction. We use a sign change at the end of the function
%       to correct for that.
%
%   See also PREPROCESS, GAUSS_CONV
%% AUTHOR    : Berook Alemayehu and Jeffrey Chiou
%% $DATE     : 04-Jun-2013 16:19:04 $
%% $Revision : 2.00 $
%% DEVELOPED : 8.0.0.783 (R2012b)
%% FILENAME  : find_fix.m

% ---------------------------BEGIN Input Parsing---------------------------
% The following 2 blocks of code enables the use of parameters, checks
% their values, and sets default values.
% You may change default values here.
p = inputParser;
p.addRequired('Data',@isstruct);

% Default: don't show any plots.
p.addParameter('showPlot',0,@(x) islogical(x)||x==0||x==1);

% Default figure save path is a new directory Figs/NAMEDATE/
p.addParameter('savePath',['Figs/' strcat(Data(1).Overview.subjectName,...
    sprintf('%1.0f',Data(1).Overview.date)) '/'],@ischar);
p.addParameter('saveFigs',1,@(x) islogical(x)||x==0||x==1);
p.addParameter('mmConvert',0.4,@(x) isnumeric(x) && isscalar(x));
p.addParameter('degConvert',330,@(x) isnumeric(x) && isscalar(x));
p.addParameter('debug',0,@(x) islogical(x)||x==0||x==1);
p.parse(Data,varargin{:});

Data = p.Results.Data;
showPlot = p.Results.showPlot;
savePath = p.Results.savePath;
saveFigs = p.Results.saveFigs;
mmConvert = p.Results.mmConvert;
degConvert = p.Results.degConvert;
debug = p.Results.debug;
% ---------------------------END Input Parsing---------------------------

% %create path to save data in
% % mkdir(savePath);

%TEMPLATE FOR SAVING FIGURES - REMOVE LATER!
% if (saveFigs==1)
%     saveas(h1,[savePath thisname '.png']);
%     saveas(h1,[savePath thisname '.fig']);
% end

% Check whether eye position has been extracted & calibrated
% Ensures preprocessing has occured successfully, particularly for
% fields we are interested in.
if ~isfield([Data(10).TrialData],'startTarget')
    Data = preprocess(Data);
end

% Initialize eyeTrials array, containing all the trials we want to analyze.
eyeTrials = [];
% DISPTEXT = true;

% Check every trial, and if the EyePosition field exists, add the trial
% index to eyeTrials. If it doesn't exist near the early trials, it is
% likely because of calibration.
for i=1:length(Data)
    if isfield([Data(i).TrialData],'EyePosition')
        if ~isempty(Data(i).TrialData.EyePosition.raw)
            eyeTrials = [eyeTrials i];
        end
    end
end

% ----------------Remaining Code is Per-Trial Analysis Code----------------
for trial = 1:length(Data)
    % This is to show how far along the code is while running.
    percentComplete(trial,length(Data)); 
   
    % For every eyeTrial
    if ismember(trial,eyeTrials)
        
        % Find raw eye position data
        rawx = Data(trial).TrialData.EyePosition.raw(:,1);
        rawy = Data(trial).TrialData.EyePosition.raw(:,2);
        
        % Convert to millimeters
        % Assuming startTarget(1) is x, (2) is y, and (3) is z.
        if str2double(Data(1).Overview.date) < 20130606
            tmpconv = 0.4;
        else
            tmpconv = 0.6;
        end
        xmm = (rawx-Data(trial).TrialData.startTarget(1))*tmpconv;
        ymm = (rawy-Data(trial).TrialData.startTarget(2))*tmpconv;
        xtargmm = (Data(trial).TrialData.endTarget(1)-Data(trial).TrialData.startTarget(1))*tmpconv;
        ytargmm = (Data(trial).TrialData.endTarget(2)-Data(trial).TrialData.startTarget(2))*tmpconv;
        
        % Working backwards, it seems 1 raw unit = 0.4 mm and 2.5 raw
        % units = 1mm
        
        % Convert to degrees. Assuming that the eye tracker is 330mm
        % away from subjects' eyes.
        xdeg = atand(xmm/330); %in degrees, where 0deg is centered on start targ
        ydeg = atand(ymm/330);
        xtargdeg = atand(xtargmm/330);
        ytargdeg = atand(ytargmm/330);
        
        
        % Review gaussian smoothing code later
        xsmooth = gauss_conv(xdeg',18); %smooth eye posn (otherwise get flickering velocity)
        ysmooth = gauss_conv(ydeg',18);
        
        if debug
            close all;
            figure; hold on;
            plot(xsmooth,'r','LineWidth',3);
            plot(xdeg);
            plot(ysmooth,'k','LineWidth',3);
            plot(ydeg,'c');
            legend('xsmooth','x','ysmooth','y');
        end
        
        % If either xsmooth or ysmooth are all NaNs, skip the rest of
        % the code
        if isnan(nanmean(xsmooth)) || isnan(nanmean(ysmooth))
            continue
        end
        
        % Find the differences to calculate velocity. Note that raw data is
        % position with ms resolution, so you get diff/ms.
        diffx = diff(xsmooth);
        diffy = diff(ysmooth);
        
        % Use sqrt of differences^2 for 2D velocity. Multiply by 1000 for
        % degrees/s
        velocity = sqrt(diffx.^2 + diffy.^2)*1000; % Actually speed, not velocity.
        
        % If velocity has no finite values (non-NaN), skip the rest of the
        % code.
        if isempty(find(~isnan(velocity)))
            continue;
        end
        
        %Threshold for determining a fixation in degrees/s
        thresh = 50;
%         thresh = 2.4*nanstd(velocity);
        
        %handle NaNs if present
        nan_ind = []; % chunks = [];
        if ~isempty(find(isnan(velocity), 1))
            
            % Find all the indices of 'velocity' where values are NaN
            nan_ind = find(isnan(velocity));
            
            % Find the finite value indices for starts and ends of NaN
            % chunks/gaps. diff(x)>1 returns a logical: 1 if x(n+1)-x(n)>1 is
            % true, 0 if false.
            
            % Starts: 'find' returns the indices of 1s in the logical noted
            % above. Add 1 to get the indices of NaNs that begin a NaN
            % chunk/gap. nan_ind(find(...)) gets the indices which correspond
            % to NaNs beginning a gap in 'velocity'. Subtract 1 for the finite
            % values. Ex:   (2  3 [NaN] NaN 4) ->
            %               (2 [3] NaN  NaN 4)
            starts = nan_ind(find(diff(nan_ind)>1)+1)-1;
            starts = [nan_ind(1)-1 starts]; % add start of 1st chunk/gap
            
            % Ends: Simpler than starts - nan_ind(...) gets the indices which
            % correspond to NaNs ending a gap in 'velocity'. Add one to get a finite
            % value. Ex:    (2 3 NaN [NaN] 4) ->
            %               (2 3 NaN  NaN [4])
            ends = nan_ind(diff(nan_ind)>1)+1;
            ends = [ends nan_ind(end)+1]; % add end of last chunk/gap
            
            % If there are NaNs at the first index of 'velocity', convert this first
            % NaN to the first gap's end. Make the first start 1 (the first
            % index).
            if(nan_ind(1)==1)
                velocity(1) = velocity(ends(1));
                starts(1) = 1;
            end
            % If there are NaNs at the last index of 'velocity', convert this first
            % NaN to the first gap's end. Make the last 'end' the last index.
            if(nan_ind(end)==length(velocity))
                velocity(end) = velocity(starts(end));
                ends(end) = length(velocity);
            end
            % The previous two blocks of code will ensure a horizontal line
            % during interpolation. For example:
            % (NaN NaN 5 6 7 NaN NaN NaN) ->
            % ( 5  NaN 5 6 7 NaN NaN  7 )
            
            % Interp1 enables both linear and other interpolation methods,
            % such as cubic and spline. This is about 6-9x slower than the
            % alternate implementation for linear interpolation, but the
            % timescale is so small it probably doesn't matter (1/2 ms)
            velocity(isnan(velocity)) = interp1(find(~isnan(velocity)),...
                velocity(~isnan(velocity)), find(isnan(velocity)),'linear');
        end
        
        % Initialize variables related to fixation.
        % thresh = 20% of highest velocity
        % min = closest local minimum (on the right side, of course)
        % ind = array of fixation start and end indices (chosen from thresh
        % and min.
        % Note: fix_start is the end of a saccade, and fix_end is the
        % beginning of a saccade.
        fix_start_thresh = []; fix_start_min = []; fix_start_ind = [];
        fix_end_thresh = []; fix_end_min = []; fix_end_ind = []; fix_time = [];
        
        % Find saccade peaks using built-in function. Finds only peaks
        % greater than 'thresh'. Note that the very first and last indices
        % of velocities will never be found as peaks (tested). Peaks have
        % to be larger than two neighboring values, and the first and last
        % only have one neighboring value.
        [peak_value, peak_ind] = findpeaks(velocity,'MINPEAKHEIGHT',thresh);
        
        % Find the threshold crossings (in ms) by taking
        % Ex. Vel = 30, thresh = 50. sign(30-50) = -1.
        % Ex2. x = diff([-1 -1 1 1 0 0 1]) = [0 2 0 -1 0 1]
        threshCrossings = diff(sign(velocity - thresh));
        % Ex3: find(x~=0) = 2 4 6
        % Note: switched to ~= from (threshCrossings>0 | threshCrossings<0)
        % for readability. Results in index, which is in ms. Bonus: faster.
        threshCrossings = find(threshCrossings~=0);
        
        %if there were NaNs in the input, put them back in so we're not creating data
        if exist('nan_ind','var')
            velocity(nan_ind) = NaN;
        end
        
        %             lostEye = find(isnan(velocity));
        %             lostEye_peaks = find(velocity(nan_ind(chunks)) > thresh);
        %             if ~isempty(lostEye_peaks)
        % %                 chunks = diff(lostEye) > 1;
        %                 lostEye_peaks = find(velocity(chunks) > thresh); %[lostEye(1)-1, find(diff(lostEye) > 1 & velocity(lostEye) > thresh), find(diff(lostEye) > 1 & velocity(lostEye) > thresh)+1 lostEye(end)+1];
        %                 if lostEye_peaks(1) == 0
        %                     lostEye_peaks = lostEye_peaks(3:end);
        % %                     lostEye = lostEye(3:end);
        %                 end
        %                 if ~isempty(lostEye_peaks) && lostEye_peaks(end) == length(velocity)
        %                     lostEye_peaks = lostEye_peaks(1:end-2);
        % %                     lostEye = lostEye(1:end-2);
        %                 end
        %                 peak_ind =  sort([peak_ind lostEye_peaks]);
        % %                 peak_ind = peak_ind(find(peak_ind));
        %                 peak_value = velocity(peak_ind);
        %             end
        
        % look at each threshold crossing
        for i = 1:length(threshCrossings)
            %-----------------Code for very first crossing-----------------
            if i == 1
                % See if there are any peaks between index 1 and the first
                % threshold crossing (index 1 can never be a peak). If
                % there are, we will need to find a fixation start point
                % (other than index/time 1).
                % Such a peak only occurs when a saccade begins before
                % recording starts, resulting in above-threshold values at
                % time 1. As the saccade is ending, it creates the first
                % threshold Xing.
                peaks = find(peak_ind > 1 & peak_ind < threshCrossings(i));
                if length(peaks) >= 1 % If there's 1 or more peaks
                    % Find the time of the first vel value from Xing to the
                    % end that is less than the 20% of the last peak value.
                    % Add the time of thresh Xing since the time found is
                    % relative to the thresh Xing.
                    fix_start_thresh = find(velocity(threshCrossings(i):end) < peak_value(peaks(end))*.20,1,'first')+threshCrossings(i);
                    % Find the time of the first local minimum after the
                    % threshold crossing
                    % Another way:
                    
                    %[temp, fix_start_min] = findpeaks(-velocity(threshCrossings(i):end));
                    % fix_start_min = fix_start_min(1)+threshCrossings(i);
                    
                    fix_start_min = find(diff(sign(diff(velocity(threshCrossings(i):end)))),1,'first')+1+threshCrossings(i);
                    
                    % Whichever one comes earlier is set as fixation start.
                    fix_start_ind = [fix_start_ind; min(fix_start_thresh,fix_start_min)];
                    % otherwise if there are no NaNs b/w beginning and the
                    % first Xing, the start will be 1. But do we need more?
                    % What if there is a below threshold set of NaNs? Do we
                    % just ignore it?
                elseif isempty(find(isnan(velocity(1:threshCrossings(i)))))
                    fix_start_ind = [fix_start_ind; 1];
                end
                
                if length(threshCrossings) == 1
                    fix_end_ind = peak_ind;
%                     peaks = peak_ind;
                else
                    % Find peaks between 1st & 2nd Xings
                    peaks = find(peak_ind > threshCrossings(i) & peak_ind < threshCrossings(i+1));
                    
                    if length(peaks) > 1        % If there are 2 or more peaks,
                        % find thresh based on velocity of first peak and min
                        if isempty(find(isnan(velocity(1:threshCrossings(i)))))
                            fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks(1))*.20,1,'last');
                            fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                            
                            if isempty(fix_end_thresh) && isempty(fix_end_min)
                                fix_end_ind = [fix_end_ind; 1];
                            elseif isempty(fix_end_min)
                                fix_end_ind = [fix_end_ind; fix_end_thresh];
                            elseif isempty(fix_end_thresh)
                                fix_end_ind = [fix_end_ind; fix_end_min];
                            else
                                fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                            end
                        end
                        
                        fix_start_thresh = find(velocity(threshCrossings(i+1):end) < peak_value(peaks(end))*.20,1,'first')+threshCrossings(i+1);
                        fix_start_min = find(diff(sign(diff(velocity(threshCrossings(i+1):end)))),1,'first')+1+threshCrossings(i+1);
                        fix_start_ind = [fix_start_ind; min(fix_start_thresh,fix_start_min)];
                        
                    elseif length(peaks) == 1       % If there is only 1 peak
                        if isempty(find(isnan(velocity(1:threshCrossings(i)))))
                            fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks)*.20,1,'last');
                            fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                            
                            if isempty(fix_end_thresh) && isempty(fix_end_min)
                                fix_end_thresh = 1;
                                fix_end_min = 1;
                            elseif isempty(fix_end_thresh)
                                fix_end_thresh = 1;
                            elseif isempty(fix_end_min)
                                fix_end_min = 1;
                            end
                            
                            fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                        end
                        
                        fix_start_thresh = find(velocity(threshCrossings(i+1):end) < peak_value(peaks)*.20,1,'first')+threshCrossings(i+1);
                        fix_start_min = find(diff(sign(diff(velocity(threshCrossings(i+1):end)))),1,'first')+1+threshCrossings(i+1);
                        if isempty(fix_start_thresh) || isempty(fix_start_min)
                            continue
                        else
                            fix_start_ind = [fix_start_ind; min(fix_start_thresh,fix_start_min)];
                        end
                    end
                end
                % -----------------Code for very last crossing-----------------
            elseif i == length(threshCrossings)
                peaks = find(peak_ind < length(velocity) & peak_ind > threshCrossings(i));
                if length(peaks) > 1
                    fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks(1))*.20,1,'last');
                    fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                    
                    if isempty(fix_end_thresh) && isempty(fix_end_min)
                        fix_end_thresh = 1;
                        fix_end_min = 1;
                    elseif isempty(fix_end_thresh)
                        fix_end_thresh = 1;
                    elseif isempty(fix_end_min)
                        fix_end_min = 1;
                    end
                    
                    fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                elseif length(peaks) == 1
                    fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks)*.20,1,'last');
                    fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                    fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                else
                    fix_end_ind = [fix_end_ind; length(velocity)];
                end
                %---------Code for crossings that aren't first or last---------
            else
                peaks = find(peak_ind > threshCrossings(i) & peak_ind < threshCrossings(i+1));
                if length(peaks) > 1
                    fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks(1))*.20,1,'last');
                    fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                    
                    if isempty(fix_end_thresh) && isempty(fix_end_min)
                        fix_end_thresh = threshCrossings(i);
                        fix_end_min = threshCrossings(i);
                    elseif isempty(fix_end_thresh)
                        fix_end_thresh = threshCrossings(i);
                    elseif isempty(fix_end_min)
                        fix_end_min = threshCrossings(i);
                    end
                    
                    fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                    
                    fix_start_thresh = find(velocity(threshCrossings(i+1):end) < peak_value(peaks(end))*.20,1,'first')+threshCrossings(i+1);
                    fix_start_min = find(diff(sign(diff(velocity(threshCrossings(i+1):end)))),1,'first')+1+threshCrossings(i+1);
                    
                    if isempty(fix_start_thresh) && isempty(fix_start_min)
                        fix_start_thresh = threshCrossings(i+1);
                        fix_start_min = threshCrossings(i+1);
                    elseif isempty(fix_start_thresh)
                        fix_start_thresh = threshCrossings(i+1);
                    elseif isempty(fix_start_min)
                        fix_start_min = threshCrossings(i+1);
                    end
                    
                    fix_start_ind = [fix_start_ind; min(fix_start_thresh,fix_start_min)];
                elseif length(peaks) == 1
                    fix_end_thresh = find(velocity(1:threshCrossings(i)) < peak_value(peaks)*.20,1,'last');
                    fix_end_min = find(diff(sign(diff(velocity(1:threshCrossings(i))))),1,'last')+1;
                    
                    if isempty(fix_end_thresh) && isempty(fix_end_min)
                        fix_end_thresh = threshCrossings(i);
                        fix_end_min = threshCrossings(i);
                    elseif isempty(fix_end_thresh)
                        fix_end_thresh = threshCrossings(i);
                    elseif isempty(fix_end_min)
                        fix_end_min = threshCrossings(i);
                    end
                    
                    if ~isempty(fix_start_ind)
                        fix_end_ind = [fix_end_ind; max(fix_end_thresh,fix_end_min)];
                    end
                    
                    fix_start_thresh = find(velocity(threshCrossings(i+1):end) < peak_value(peaks)*.20,1,'first')+threshCrossings(i+1);
                    fix_start_min = find(diff(sign(diff(velocity(threshCrossings(i+1):end)))),1,'first')+1+threshCrossings(i+1);
                    if isempty(fix_start_thresh)
                        fix_start_thresh = threshCrossings(i+1);
                    end
                    if isempty(fix_start_min)
                        fix_start_min = threshCrossings(i+1);
                    end
                    fix_start_ind = [fix_start_ind; min(fix_start_thresh,fix_start_min)];
                end
            end
        end
        
        %        removeFixation = find([fix_end_ind-fix_start_ind] < 0);
        %             fix_start_ind(removeFixation) = [];
        %             fix_end_ind(removeFixation) = [];
        
        
        %         if trial == 500
%         subplot(2,1,2); plot(velocity);
%         xlabel('Time(ms)');
%         ylabel('Eye Velocity (degrees/s)');
%         hold on
%         if Data(trial).Overview.trialStatus == 1
%             line([Data(trial).TrialData.timeDelayStart_host Data(trial).TrialData.timeDelayStart_host],ylim,'Color','k','LineStyle','--');
%             line([Data(trial).TrialData.timeGoCue_host Data(trial).TrialData.timeGoCue_host],ylim,'Color','k','LineStyle','--');
%         end
%         for g = 1:length(fix_start_ind)
%             subplot(2,1,2); hold on; line([fix_start_ind(g) fix_start_ind(g)],ylim,'Color','g','LineStyle','--');
%             subplot(2,1,2); hold on; line([fix_end_ind(g) fix_end_ind(g)],ylim,'Color','r','LineStyle','--');
%         end
%         subplot(2,1,1); plot(xsmooth); hold on; plot(ysmooth,'r');
%         ylabel('Eye Position (degrees)');
%         for g = 1:length(fix_start_ind)
%             subplot(2,1,1); hold on; line([fix_start_ind(g) fix_start_ind(g)],ylim,'Color','g','LineStyle','--');
%             subplot(2,1,1); hold on; line([fix_end_ind(g) fix_end_ind(g)],ylim,'Color','r','LineStyle','--');
%         end
%         subplot(2,1,2); hold on; line(xlim,[thresh thresh],'Color','k','LineStyle','--');
%         subplot(2,1,2); hold on; line(xlim,[nanstd(velocity)*2.4 nanstd(velocity)*2.4],'Color','k','LineStyle','--');
%         saveas(gca,['S:\Animals\Ford\Analysis Results\Fixation Behavior\Eye Traces\eyeTrace_Trial ' num2str(trial)],'pdf')
%         close all
        %         end
        
        
        if showPlot
            % close all;
            figure;
            subplot(2,1,1)
            plot(xdeg)
            hold on;
            plot(ydeg,'g')
            ylabel('Degrees')
%             ylim = get(gca,'YLim');
            line([Data(trial).TrialData.timeTargetOn_photo Data(trial).TrialData.timeTargetOn_photo], ylim,'Color', [0 0 0])
            line([Data(trial).TrialData.timeTargetAcquired_photo Data(trial).TrialData.timeTargetAcquired_photo], ylim,'Color', [0 0 0])
            title(trial)
            
            subplot(2,1,2)
            plot(velocity,'k')
            hold on
            line([0 length(velocity)], [thresh,thresh])
            ylabel('Degrees/s')
            %plot(locs,pks,'r*')
        end
        
        
        if ~isempty(fix_end_ind)
            xposdeg = []; yposdeg = []; xposmm = []; yposmm = []; fix_time = [];
            for i=1:length(fix_start_ind)
                this_fix = fix_start_ind(i):fix_end_ind(i);
                
                if strcmp(Data(1).Overview.subjectName,'Ford')
                    %fixation has to be at least 100 ms long
                    if length(this_fix) >= 100
                        this_xposdeg = nanmean(xdeg(this_fix));
                        this_yposdeg = nanmean(ydeg(this_fix));
                        this_xposmm = nanmean(xmm(this_fix));
                        this_yposmm = nanmean(ymm(this_fix));
                        
                        xposdeg = [xposdeg; this_xposdeg];
                        yposdeg = [yposdeg; this_yposdeg];
                        xposmm = [xposmm; this_xposmm];
                        yposmm = [yposmm; this_yposmm];
                        fix_time = [fix_time; [this_fix(1), this_fix(end)]];
                        
                        if showPlot
                            line([this_fix(1), this_fix(1)],[0 max(velocity)],'Linestyle','--','Color',[1 0 0])
                            line([this_fix(end), this_fix(end)],[0 max(velocity)],'Linestyle','--','Color',[0 1 0])
                        end
                    end
                else
                    if length(this_fix) >= 150
                        this_xposdeg = nanmean(xdeg(this_fix));
                        this_yposdeg = nanmean(ydeg(this_fix));
                        this_xposmm = nanmean(xmm(this_fix));
                        this_yposmm = nanmean(ymm(this_fix));
                        
                        xposdeg = [xposdeg; this_xposdeg];
                        yposdeg = [yposdeg; this_yposdeg];
                        xposmm = [xposmm; this_xposmm];
                        yposmm = [yposmm; this_yposmm];
                        fix_time = [fix_time; [this_fix(1), this_fix(end)]];
                        
                        if showPlot
                            line([this_fix(1), this_fix(1)],[0 max(velocity)],'Linestyle','--','Color',[1 0 0])
                            line([this_fix(end), this_fix(end)],[0 max(velocity)],'Linestyle','--','Color',[0 1 0])
                        end
                    end
                end
            end
            
            %save to the data struct for later analysis
            Data(trial).TrialData.EyePosition.xfixdeg = xposdeg; %in degrees
            Data(trial).TrialData.EyePosition.yfixdeg = -yposdeg; %y is neg b/c host graphs coords upside down
            Data(trial).TrialData.EyePosition.xfixmm = xposmm;
            Data(trial).TrialData.EyePosition.yfixmm = -yposmm;
            Data(trial).TrialData.EyePosition.fixTimes = fix_time; %time from trial start
            Data(trial).TrialData.EyePosition.xtargmm = xtargmm;
            Data(trial).TrialData.EyePosition.ytargmm = -ytargmm;
            Data(trial).TrialData.EyePosition.xtargdeg = xtargdeg;
            Data(trial).TrialData.EyePosition.ytargdeg = -ytargdeg;
        end
        
        %         if isempty(xposmm)
        %             keyboard;
        %         end
        
        if isempty(fix_time) || isempty(fix_end_ind)
            fprintf(['Error: No Valid Fixations for Trial ', num2str(trial) '\n']);
            Data(trial).TrialData.EyePosition.xfixdeg = NaN; %in degrees
            Data(trial).TrialData.EyePosition.yfixdeg = NaN;
            Data(trial).TrialData.EyePosition.xfixmm = NaN; %in mm
            Data(trial).TrialData.EyePosition.yfixmm = NaN;
            Data(trial).TrialData.EyePosition.fixTimes = [NaN NaN];
            %             Data(trial).TrialData.EyePosition.xtargmm = NaN;
            %             Data(trial).TrialData.EyePosition.ytargmm = NaN;
            %             Data(trial).TrialData.EyePosition.xtargdeg = NaN;
            %             Data(trial).TrialData.EyePosition.ytargdeg = NaN;
        end
    end
    %     end
end

end