clear; clc; clf; close all

%Create signal
figure
binWidth = 1; %ms
numSamp = 1000;
x = ones(1,numSamp);
x = sin(1:numSamp);
x(100:150) = -1;
x(700:800) = 2;
plot(x)

%Create gaussian kernel with sigma = 25ms
sigma = 25/binWidth; %sigma in number of bins 
windowSideLength = ceil(3*sigma); %half width of gaussian kernal in bins
kernel = normpdf([-windowSideLength:1:windowSideLength],0,sigma);
kernel = kernel.*(1/sum(kernel));
figure
plot(kernel)

%Convolve
c = nan(1,length(x));
kernelMidInd = median(1:length(kernel));
kernelLength = length(kernel);

for i = windowSideLength+1:length(x)-windowSideLength
    tempSum = 0;
    for j = i-windowSideLength:i+windowSideLength
       tempSum = x(j)*kernel(kernelLength-(j-(i-windowSideLength))) + tempSum; 
%        j
%        i-j+windowSideLength+1
    end
    c(i) = tempSum;
end


%%Convolve w matlab code 
cMat = nan(1,length(x));
cMat(1,windowSideLength+1:length(x)-windowSideLength) = conv(x,kernel,'valid');

figure
plot(c)
hold on
plot(cMat)
% for sort = 1:size(X,2)
%     w1 = conv(X(:,sort),kernel,'same');
%     w2 = conv(X(:,sort),kernel);
%     StatesTraj =1
% end