%% Calculates R^2 value, given Y (1D target variable) and Ypred (model predictions for Y)

function [r2] = getR2(Y,Ypred)

    N = length(Y);
    
    %Get SSE
    SSE = 0;
    for i = 1:N
       SSE = (Y(i)-Ypred(i)).^2 + SSE; 
    end
    %Get SST
    SST = 0;
    u = mean(Y);
    for i = 1:N
        SST = (Y(i)-u).^2 + SST;
    end
    
    %Get R^2
    r2 = 1 - SSE/SST;

end