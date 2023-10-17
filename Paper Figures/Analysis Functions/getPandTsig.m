function [pSig,tSig] = getPandTsig(trajStruct,field,minNumTimestamps,postureList,numPostures,targetList,numTargets)
        %Get posture and target components of neural activity. Field
        %specifies whether to use z-scored, smoothed FR or PCA scores

        %Form X
        X = getX(trajStruct,field,minNumTimestamps,postureList,numPostures,targetList,numTargets);    

        % Do marginalizations of X (Xdims: 1=time, 2=target, 3=posture, 4=channel)
        %Condition-invariant
        CIMargOffset = mean(X,[1 2 3],'omitnan');
        CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
        %Posture and Target Traj
        tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
        pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

end