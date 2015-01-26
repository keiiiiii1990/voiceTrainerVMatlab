function f0Structure = extractInitialF0byHigherSymmetry(x,fs,opt)
%   f0Structure = extractInitialF0byHigherSymmetry(x,fs)

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
if nargin == 3
    if isfield(opt,'f0floor')
        f0floor = opt.f0floor;
    else
        opt.f0floor = 14; % revision for extracting 32Hz f0
    end;
    if isfield(opt,'f0ceil')
        f0ceil = opt.f0ceil;end;
    if isfield(opt,'nInOctave')
        nInOctave = opt.nInOctave;end;
    if isfield(opt,'frameShift')
        frameShift = opt.frameShift;end;
else
    opt.f0floor = 14; % revision for extracting 32Hz f0
end;

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