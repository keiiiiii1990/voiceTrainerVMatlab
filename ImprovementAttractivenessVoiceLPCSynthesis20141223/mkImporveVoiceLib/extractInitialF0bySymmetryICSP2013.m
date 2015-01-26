function f0Structure = extractInitialF0bySymmetryICSP2013(x,fs)
%   f0Structure = extractInitialF0bySymmetryICSP2013(x,fs)

%   02/Sep./2012 by Hideki Kawahara
%   13/Sep./2012 Optimized for ICASSP 2013 architecture
%   23/Sep./2012 f0floor defect fixed

%%
startTime = tic;
w = [1 64 128];
w = w/sum(w);
magnifier = 4;
exponent = 4;
normalizer = [0.3472 0.2308 0.0720];
opt.f0floor = 14; % revision for extracting 32Hz f0

r = f0bySymmetryV4(x,fs,opt);
r.updatedPeriodicityMap = ...
    individualToMixedMeasure(r,normalizer,w,magnifier,exponent);
candidates = candidatesFromTwoMaps(r);

f0Structure = r;
f0Structure.candidates = candidates;
f0Structure.f0Raw = candidates.f0Initial(1,:)';
f0Structure.rawPeriodicity = candidates.periodicityList(1,:)';
f0Structure.totalElapsedTime = toc(startTime);
%%
return;