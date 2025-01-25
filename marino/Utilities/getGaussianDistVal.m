function [p] = getGaussianDistVal(x,u,S)
x=x';
u=u';
 p = (det(2*pi*S).^-0.5)*exp(-0.5*(x-u)'*inv(S)*(x-u));


end