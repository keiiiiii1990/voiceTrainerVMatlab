function f0Struct = higherSymKalmanInitial(x,fs)
%   f0Struct = higherSymKalmanInitial(x,fs)
%

%   Designed and coded by Hideki Kawahara
%   18/Aug./2013
%   18/Nov./2013 frontend only

startTic = tic;
if size(x,1) < size(x,2);x=x';end;
x = x(:,1);
f0RawStr = extractInitialF0bySymmetryICSP2013(x,fs);
spex2 = periodicityIndex2VbasedWeight2(f0RawStr);
avF0cent = 1200*log2(spex2.averagedF0/440);
sigmaObs = spex2.totalSigma;
phi2 = 1600;
kalmanStr = kalmanF0Smoother(avF0cent,sigmaObs,phi2);
latentF0 = 440*2.0.^(kalmanStr.latentF0Cent/1200);
f0Struct.f0InitialByHigherSymmetry = f0RawStr;
f0Struct.estimatedVariances = spex2;
f0Struct.kalmanSmootherStructure = kalmanStr;
f0Struct.latentF0 = latentF0(:);
f0Struct.latentSDcoeff = 2.0.^(sqrt(kalmanStr.latentVarCent2)/1200);
f0Struct.averagedF0 = spex2.averagedF0;
f0Struct.totalObservationVariance = spex2.totalSigma;
f0Struct.temporalPositions = f0RawStr.temporalPositions;
f0Struct.elapsedTime = toc(startTic);
return;