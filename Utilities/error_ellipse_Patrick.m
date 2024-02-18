function [X,Y,Z] = error_ellipse_Patrick(data,pct)

%Generates output that allows you to draw error ellipse for data (assuming Gaussian distribution. 
%Takes as input data (must be 2 or 3d) and pct, which specifies the
%percentage of data to be contained w/in the ellipse.
%Outputs X,Y,Z.  For 2d data, X and Z allow you to draw error
%ellipse with plot(X,Y).  For 3d, use surf(X,Y,Z);

%Get dimensionality  and mean of data
d = size(data,2);
mu = mean(data);
if ~ismember(d,[2,3])
   error('data must be 2 or 3 dimensional in current implementation of function'); 
end

%Get eignenvectors and values of covariance matrix; sort
[eigVec,eigVal] =  eig(cov(data));
eigVal = diag(eigVal); 
[eigVal, sortOrder] = sort(eigVal, 'descend'); 
eigVec = eigVec(:, sortOrder);

%Get radii of ellipse axes
c = chi2inv(pct/100,d);
radii = sqrt(c*eigVal);

%If 2d - output X,Y to plot(X,Y) to draw ellipse
if d == 2
    thetaList = linspace(0,2*pi);
    G = radii(1)*cos(thetaList);
    H = radii(2)*sin(thetaList);
    rot = transpose(eigVec*vertcat(G,H));
    X = rot(:,1) + mu(1);
    Y = rot(:,2) + mu(2);
    Z = [];
end

%If 3d - output X,Y,Z for surf function to use as inputs
if d == 3
    [G,H,K] = ellipsoid_Patrick(0,0,0,radii(1),radii(2),radii(3),15);
    numRow = size(G,1); numCol = size(G,2);
    X = zeros(numRow,numCol); Y = zeros(numRow,numCol); Z = zeros(numRow,numCol);
    for i = 1:numRow
        for j = 1:numCol
            X(i,j) = eigVec(1,:)*[G(i,j);H(i,j);K(i,j)] + mu(1);
            Y(i,j) = eigVec(2,:)*[G(i,j);H(i,j);K(i,j)] + mu(2);
            Z(i,j) = eigVec(3,:)*[G(i,j);H(i,j);K(i,j)] + mu(3);
        end
    end
end 

end