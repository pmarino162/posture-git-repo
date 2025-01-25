%Performs F test on 
function [F,p] = myFTest(Y,YpredF)

n = size(Y,1);
YpredR = ones(n,1)*mean(Y);
SSER = 0; SSEF = 0;
for i = 1:n
    SSER = (Y(i,1)-YpredR(i,1)).^2 + SSER;

    SSEF = (Y(i,1)-YpredF(i,1)).^2 + SSEF;
end

dfR = n-1;

dfF = n-3;

F = ((SSER-SSEF)/(dfR-dfF))/(SSEF/dfF);

p = fpdf(F,dfR-dfF,dfF);

end