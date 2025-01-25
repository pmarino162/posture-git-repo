clear; clc; close all
day = [1:13,18 20 21];
chOff = [0 0 0 4 4 6 10 10 12 21 19 16 14 15 18 15];
M1ChOff = [0 0 0 3 3 4 8 9 11 19 17 15 13 12 13 12];

plot(day,chOff,'k')
hold on
plot(day,chOff,'.','MarkerSize',12,'Color','k')

plot(day,M1ChOff,'b')
hold on
plot(day,M1ChOff,'.','MarkerSize',12,'Color','b')

ylabel('# Channels Off')
xlabel('Day')
