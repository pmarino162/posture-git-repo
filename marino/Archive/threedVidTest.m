clear; clc; close all

data = normrnd(0,5,100,3);
figure
plot3(data(:,1),data(:,2),data(:,3),'.','MarkerSize',15);
xlabel('x'); ylabel('y');zlabel('z');

OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;

CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)
