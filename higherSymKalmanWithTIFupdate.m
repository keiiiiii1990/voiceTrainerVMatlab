function f0Struct = higherSymKalmanWithTIFupdate(x,fs)
%   f0Struct = higherSymKalmanWithTIFupdate(x,fs)
%

%   Designed and coded by Hideki Kawahara
%   18/Aug./2013

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
periodicityStructure.temporalPositions = f0RawStr.temporalPositions;
periodicityStructure.f0Raw = latentF0(:);
refinedStr = refineExGxOutputRevSqx(x,fs,periodicityStructure);%,optStr);
f0Struct.f0InitialByHigherSymmetry = f0RawStr;
f0Struct.estimatedVariances = spex2;
f0Struct.kalmanSmootherStructure = kalmanStr;
f0Struct.refinedF0byTIFStructure = refinedStr;
f0Struct.F0 = refinedStr.f0Refined;
f0Struct.latentF0 = latentF0;
f0Struct.latentSDcoeff = 2.0.^(sqrt(kalmanStr.latentVarCent2)/1200);
f0Struct.averagedF0 = spex2.averagedF0;
f0Struct.totalObservationVariance = spex2.totalSigma;
f0Struct.temporalPositions = refinedStr.temporalPositions;
f0Struct.elapsedTime = toc(startTic);
return;