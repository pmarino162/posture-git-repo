% Gets parameters for trajStruct function used in most analyses of posture paper

function [condFields,trajFields,trialInclStates,binWidth,kernelStdDev] = getTrajStructParamsKinematicComparison(dataset)

        binWidth = 25;
        kernelStdDev = 25;
        
        
        trialInclStates = struct('trialName','','inclStates',[]);
        
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                trajFields = {'zSmoothFR','bciCursorTraj'};
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
                task = 'bci';
            case {'N20171215','N20180221'}
                trajFields = {'zSmoothFR','bciCursorTraj'};
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Release','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trajFields = {'zSmoothFR','bciCursorTraj'};
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';

            %Iso
            case {'E20200116','E20200117','E20200120'}
                trajFields = {'zSmoothFR','force'};
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
                task = 'iso';
                
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trajFields = {'zSmoothFR','markerPos'};
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target Acquire','first',0},{'state','Target Hold','first',0}};
                task = 'reach';
                
            case {'R20200221','R20200222'}
                trajFields = {'zSmoothFR','cursorPos'};
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
                task = 'reach';
                        
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trajFields = {'zSmoothFR','cursorPos'};
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
                task = 'reach';
        end 

end