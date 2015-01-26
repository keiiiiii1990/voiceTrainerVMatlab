function kalmanStr = kalmanF0Smoother(avF0cent,sigmaObs,phi2)
%   kalmanStr = kalmanF0Smoother(avF0cent,sigmaObs,phi2)

%   Designed and coded by Hideki Kawahara
%   15/Aug./2013 

startTic = tic;
nFrames = length(avF0cent);

sigma02 = 2400^2;
mu0 = 0;
Mplus = zeros(nFrames,1);
Vplus = zeros(nFrames,1);
Pprd = zeros(nFrames+1,1);

Mplus(1) = (avF0cent(1)*sigma02+mu0*sigmaObs(1))/(sigma02+sigmaObs(1));
Vplus(1) = (sigma02*sigmaObs(1))/(sigma02+sigmaObs(1));
Pprd(2) = phi2+Vplus(1);

for ii = 2:nFrames
    Mplus(ii) = (avF0cent(ii)*Pprd(ii)+Mplus(ii-1)*sigmaObs(ii))/(Pprd(ii)+sigmaObs(ii));
    Vplus(ii) = (Pprd(ii)*sigmaObs(ii))/(Pprd(ii)+sigmaObs(ii));
    Pprd(ii+1) = phi2+Vplus(ii);
end;

Mminus = zeros(nFrames,1);
Vminus = zeros(nFrames,1);

Mminus(nFrames) = Mplus(nFrames);
Vminus(nFrames) = Vplus(nFrames);

for ii = nFrames:-1:2
    Mminus(ii-1) = (Vplus(ii-1)*Mminus(ii))/(phi2+Vplus(ii-1))+(Mplus(ii-1)*phi2)/(phi2+Vplus(ii-1));
    Vminus(ii-1) = Vplus(ii-1)/(phi2+Vplus(ii-1))*(phi2+(Vplus(ii-1)*Vminus(ii))/(phi2+Vplus(ii-1)));
end;
kalmanStr.elapsedTime = toc(startTic);
kalmanStr.latentF0Cent = Mminus;
kalmanStr.latentVarCent2 = Vminus;
return;