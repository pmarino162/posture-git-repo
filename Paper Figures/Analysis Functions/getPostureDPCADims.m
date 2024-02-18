function [W,V,postureDims,targetDims,explVar,additionalVarExpl] = getPostureDPCADims(trajStruct,regularize,postureList,numPostures,targetList,numTargets,numChannels,minNumTimestamps)

        %Applies dPCA and returns the posture and target dims
        %used for visualizations in the paper. explVar is ordered: posture
        %1, posture 2, target 1, target 2

        %% Form X and permute to match dPCA paper
        [X] = getX(trajStruct,'avgZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets);    
        Xdpca = permute(X,[4,3,2,1]);
                
        %% dPCA Parameters  
        numDPCs = 15;
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        N = numChannels; P = numPostures; T = numTargets;
        time = trajStruct(1).avgZSmoothFR.timestamps(1:minNumTimestamps);
        timeEvents = [];

        %% Run dPCA
        %Regularalized version
        if regularize
            %Get Xfull and permute
            [Xfull] = getXFull(trajStruct,'allZSmoothFR',minNumTimestamps,postureList,numPostures,targetList,numTargets); 
            XdpcaFull = permute(Xfull,[4,3,2,1,5]);
            %Set up 'numOfTrials' for use with dpca_optimizeLambda
            numOfTrials = zeros(size(Xdpca,1),size(Xdpca,2),size(Xdpca,3));
            for i = 1:size(Xdpca,1)
                for j = 1:size(Xdpca,2)
                    for k = 1:size(Xdpca,3)
                        numOfTrials(i,j,k) = sum(~isnan(XdpcaFull(i,j,k,1,:)));
                    end
                end
            end
            %Compute optimal lambda
            optimalLambda = dpca_optimizeLambda(Xdpca, XdpcaFull, numOfTrials, ...
                'combinedParams', combinedParams, ...
                'simultaneous', true, ...
                'numRep', 10, ...  % increase this number to ~10 for better accuracy
                'filename', 'tmp_optimalLambdas.mat');
            %Compute Cnoise
            Cnoise = dpca_getNoiseCovariance(Xdpca, ...
                XdpcaFull, numOfTrials, 'simultaneous', true);
            % Compute dPCs
            [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
                'combinedParams', combinedParams, ...
                'lambda', optimalLambda, ...
                'Cnoise', Cnoise);
        
        %Unregularized version
        else
            [W,V,whichMarg] = dpca(Xdpca, numDPCs, ...
                'combinedParams', combinedParams);
        end
        
        %% Get variance explained by dPCs
        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);
        cumulativeDPCA = explVar.cumulativeDPCA;
        additionalVarExpl = [cumulativeDPCA(1) diff(cumulativeDPCA)];

        %% Get IDs of dPCs pertaining to posture, target, and CI (unused)
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);

        
end