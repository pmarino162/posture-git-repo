function [PD,MD,b0,p] = fitTC(targetMeans) 

%Fits a tuning curve to 8-target CO means for 1 channel
    targetList = 1:8; numTargets = length(targetList);
    targetAngles = transpose(45*(targetList-1)); 
    y = targetMeans'; x = nan(numTargets,3);
    for targetInd = targetList
        x(targetInd,:) = [1,sind(targetAngles(targetInd)),cosd(targetAngles(targetInd))];
    end
    [B,~,~,~,stats] = regress(y,x);
    b0 = B(1,1); b1 = B(2,1); b2 = B(3,1); 
    p = stats(3); %P value from F-test
    MD = sqrt(b1.^2 + b2.^2);
   
    PD = atan2d(b1,b2);            
    if PD < 0 %Map [-180,180] to [0,360]
        PD = 360 - abs(PD);
    end
    

end