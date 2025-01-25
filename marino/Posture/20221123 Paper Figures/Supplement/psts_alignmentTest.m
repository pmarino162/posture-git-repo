clear; clc;

A = [0,0; 1,1; 2,1.5];

B = [3,3; 3,4; 5,7];

[PA,~,~,~,explA] = pca(A);

[PB,~,~,~,explB] = pca(B);

figure
scatter(A(:,1),A(:,2),10,'r')
hold on 
scatter(B(:,1),B(:,2),10,'b')

varA = trace(cov(A));
varAB = trace(cov(A*PB(:,1)));
varAA = trace(cov(A*PA(:,1)));
varBA = trace(cov(B*PA(:,1)));
varBB = trace(cov(B*PB(:,1)));

varAB/varAA
varBA/varBB