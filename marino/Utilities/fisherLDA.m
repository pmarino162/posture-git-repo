% Function to perform Fisher LDA 
%
% Usage: dirs = fisherLDA(pts, labels)
%
% Input: 
%
%   pts - a matrix of data points. Each data point is a row
%
%   labels - a vector labels for the data in pts
%
% Outputs:
%
%   dirs - a matrix of projection directions.  Each column is a direction.
%   Will return min(C - 1, D) vectors, where C is the number of classes in 
%   the data and D is the dimensionality of the data. 
%
% Author: wbishop@cs.cmu.edu
% Edited: 6/5/2020 Patrick Marino pmarino162@gmail.com

function [dirs,eigSpec] = fisherLDA(pts, labels)

dims = size(pts,2); 
uniqueLabels = unique(labels); 
nLabels = length(uniqueLabels); 

withinScatter = zeros(dims); 

classMns = nan(nLabels, dims); 
for lI = 1:nLabels
    curLabel = uniqueLabels(lI);
    classMns(lI,:) = mean(pts(labels == curLabel,:),1);
    
    withinScatter = withinScatter + ...
        cov(pts(labels == curLabel,:)); 
end
betweenScatter = cov(classMns); 

[eigVcs, eigVls] = eig(betweenScatter, withinScatter); 

eigVls = diag(eigVls); 
[eigVls, sortOrder] = sort(eigVls, 'descend'); 

eigVcs = eigVcs(:, sortOrder);

nReturnDims = min(nLabels-1, dims); 

dirs = eigVcs(:, 1:nReturnDims); 
eigSpec = eigVls(1:nReturnDims,1);
