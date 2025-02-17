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
%             case {'E20210901'} %Earl elbow and shoulder session
%                 trialInclStates(1).trialName = {'BCI Center Out'};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%                 trialInclStates(1).inclStates = {{'state','Step 1','first',50},{'state','Step 1','first',250}};
%             case {'E20190729'} %Earl alternate decoder session                
%                 trialInclStates(1).trialName = {'CentetOut_BC_TouchBar'};
%                 %trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
%                 trialInclStates(1).inclStates = {{'state','Step 1','first',50},{'state','Step 1','first',250}};
%                 trialInclStates(2).trialName = {'CenterOutCenter_BC_TouchBar'};
%                 %trialInclStates(2).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
%                 trialInclStates(2).inclStates = {{'state','Step 1','first',50},{'state','Step 1','first',250}};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'},{'decoder','conditionData','decoderPostureID'}};
%             case {'E20190830'}
%                 trialInclStates(1).trialName = {'CenterOut_ForceBar_BC'};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%                 trialInclStates(1).inclStates = {{'state','BC Freeze','first',0},{'state','Success with Reward','first',0}};
% 
%             %Iso
%             case {'E20200116','E20200117','E20200120'}
%                 trialInclStates(1).trialName = {'IsometricForce_1D'};   
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%                 trialInclStates(1).inclStates = {{'state','Target','first',50},{'state','Target','first',250}};
%                 %trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-250},{'kin','moveOnsetTime','first',0}};
%                 task = 'iso';
%             %Reach
%             case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
%                 trialInclStates(1).trialName = {'GridReaching'};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%                 trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
%                 task = 'reach';
            case {'R20200221','R20200222'}
                trajFields = {'zSmoothFR','markerPos'};
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
                task = 'reach';
%                 
%             case {'R20200221_all_visual_data','R20200222_all_visual_data'} %Includes both visual target locations
%                 trialInclStates(1).trialName = {'Rocky Dissociation'};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'},{'visual','conditionData','visualID'}};
%                 trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
%                 task = 'reach';                
%                 
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trajFields = {'zSmoothFR','markerPos'};
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Reach','first',0},{'state','Target Hold','first',0}};
                task = 'reach';
%             
%             case {'N20190222_all_visual_data','N20190226_all_visual_data','N20190227_all_visual_data','N20190228_all_visual_data','N20190307_all_visual_data'}
%                 % Includes both visual target locations
%                 trialInclStates(1).trialName = {'Nigel Dissociation'};
%                 condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'},{'visual','conditionData','visualID'}};
%                 trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
%                 task = 'reach';    
%                 
%             %Multiple tasks paradigm
%             case {'E20200314'}         
%                 condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%                 trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
%                     trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',50},{'state','Step 1 Freeze','first',250}};
%                 trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
%                     trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
%                 trialInclStates(3).trialName = {'IsometricForce_1D'};
%                     trialInclStates(3).inclStates = {{'state','Target','first',50},{'state','Target','first',250}};
%                     %trialInclStates(3).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        end 


end