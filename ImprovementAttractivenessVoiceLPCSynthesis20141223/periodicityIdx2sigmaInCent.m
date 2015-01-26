function sigmaInCent = periodicityIdx2sigmaInCent(periodicityIdx)
sigmaInCent = 2157^2.0./(1+exp(15.3171*(periodicityIdx-0.16))) ...
    +(1-1.0./(1+exp(15.3171*(periodicityIdx-0.16)))).*(120*(1-periodicityIdx)).^2;
return;