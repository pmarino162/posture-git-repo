function [numCondTrials] = getNumCondTrials(trajStruct,varargin)    
    
        %Variable arguments
        showHist = false;
        assignopts(who,varargin);
        
        %Get number of trials in each condition
        numCondTrials = nan(1,size(trajStruct,2));
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allZSmoothFR,2);
           numCondTrials(i) = numTraj;
        end
        
        %Optional histogram
        if showHist
            figure
                histogram(numCondTrials)
                xlabel('Number of trials')
                ylabel('Number of conditions')
        end
end