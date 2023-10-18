%This is copied from 'causalInferencePrecomp' in the GPFA package
clear;clc;
%Load bci params from one session
load('D:\Animals\Earl\2020\03\20200316\MkHstDec_GPFA_Earl20200316_Intuitive_15.mat')
CRinv = bci_params.CRinv;
Tstar = bci_params.Tmax; %maximum time = 22*45~=1000ms
params.C = zeros(10,10);
params.eps = zeros(1,10);
params.a = zeros(1,10);
params.gamma = 1;
params.covType = 'logexp';
for T = 1:Tstar
    [K_big, K_big_inv] = make_K_big(params, T);        
    K_big = sparse(K_big);

    blah        = cell(1, T);
    [blah{:}]   = deal(CRinvC);
    invM        = invPerSymm(K_big_inv + blkdiag(blah{:}), xDim,... 
			     'offDiagSparse', true);
    
    % Compute blkProd = CRinvC_big * invM efficiently
    % blkProd is block persymmetric, so just compute top half
    Thalf   = ceil(T/2);
    blkProd = zeros(xDim*Thalf, xDim*T);
    idx     = 1: xDim : (xDim*Thalf + 1);
    for t = 1:Thalf
      bIdx            = idx(t):idx(t+1)-1;
      blkProd(bIdx,:) = CRinvC * invM(bIdx,:);
    end
    % Compute just the first block row of blkProd
    blkProd = K_big(1:xDim, :) *... 
                fillPerSymm(speye(xDim*Thalf, xDim*T) - blkProd, xDim, T);

    % Last block row of blkProf is just block-reversal of first block row
    horzIdx = bsxfun(@plus, (1:xDim)', ((T-1):-1:0)*xDim);
    M       = blkProd(:, horzIdx(:));

    precomp.filt(T).M = M;
  end