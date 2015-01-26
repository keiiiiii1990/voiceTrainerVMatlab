%%  F0 extractor based on an extended zero crossing
function f0Structure = f0bySymmetryV4(xk,fs,opt)
%   f0Structure = f0bySymmetryV4(xk,fs,opt)

%   Old version was f0byExtendedZerocrossingFixpUDRev4(xk,fs,opt)
%   Interval based F0 extractor with symmetry measure detection
%   Designed and coded by Hideki Kawahara
%   10/Oct./2011 for extreme based
%   12/Oct./2011 debug
%   14/Oct./2011 debug
%   16/Oct./2011 otpimization based on profiling
%   07/Apr./2012 revised for Interspeech 2012
%   27/Aug./2012 higher order asymmetry
%   02/Sep./2012 re-design

%tic;
f0floor = 32;
f0ceil = 1300;%1000; % for showa
f0ceil = 1000;
nInOctave = 6;
frameShift = 0.001;

xk = xk(:);
%xkOrg = xk;
xk(xk==0) = randn(length(xk(xk==0)),1)*std(xk)*0.00001; % was serious bug. fixed

%%
if nargin == 3
    if isfield(opt,'f0floor')
        f0floor = opt.f0floor;end;
    if isfield(opt,'f0ceil')
        f0ceil = opt.f0ceil;end;
    if isfield(opt,'nInOctave')
        nInOctave = opt.nInOctave;end;
    if isfield(opt,'frameShift')
        frameShift = opt.frameShift;end;
end;
%%
tic
analysisConditions.f0floor = f0floor;
analysisConditions.f0ceil = f0ceil;
analysisConditions.nInOctave = nInOctave;
analysisConditions.frameShift = frameShift;

targetSamplingFrequency = f0ceil*8;
if floor(fs/targetSamplingFrequency) > 0
    newFs = fs/floor(fs/targetSamplingFrequency);
    %xStart = decimate(xk,floor(fs/targetSamplingFrequency)); replaced!
    newFc = newFs/2*0.9;
    transisitionWidth = newFs/2*0.1;
    windowLengthHalf = round(fs/transisitionWidth/2*4);
    lpfForDownSamplingBase = nuttallwin12(windowLengthHalf*2+1); % 2012/03/11
    timeForWindow = (-windowLengthHalf:windowLengthHalf)'/fs;
    sincFunctionForCarrier = sin(2*pi*newFc*timeForWindow+0.000000000001)./(2*pi*newFc*timeForWindow+0.000000000001);
    lpfForDownSampling = lpfForDownSamplingBase.*sincFunctionForCarrier;
    xStart = fftfilt(lpfForDownSampling,[xk;zeros(3*windowLengthHalf,1)]);
    xStart = xStart(windowLengthHalf+(1:floor(fs/targetSamplingFrequency):length(xk)));
else
    newFs = fs;
    xStart = xk;
end;
fcList = f0floor*2.0.^(0:1/nInOctave:ceil(log2(f0ceil/f0floor)*nInOctave)/nInOctave);
numberOfChannels = length(fcList);
%maxFixedP = round(numberOfChannels/2);
bestLength = round(1.5*fs/fcList(1)/2)*2+1;
ww = nuttallwin12(bestLength);
xTest = [xStart;zeros(length(ww*1.5),1)];
baseIndex = 1:length(xStart);
smoothedX = fftfilt(ww,xTest);
crossingStructure = extendedSymmetryMeasure(smoothedX((bestLength-1)/2+baseIndex),newFs,frameShift);%std(fftfilt(ww,x));
temporalPositions = crossingStructure.temporalPositions;
numberOfLocations = length(temporalPositions);
f0CandidatesMap = zeros(numberOfChannels,numberOfLocations);
smoothedAMMap = zeros(numberOfChannels,numberOfLocations);
smoothedFMMap = zeros(numberOfChannels,numberOfLocations);
smoothedSMMap = zeros(numberOfChannels,numberOfLocations);
amplitudeMap = zeros(numberOfChannels,numberOfLocations);

for ii = 1:numberOfChannels
    bestLength = round(1.5*newFs/fcList(ii)/2)*2+1;
    ww = nuttallwin12(bestLength);
    smoothedX = fftfilt(ww,xTest);
    crossingStructure = ...
        extendedSymmetryMeasure(smoothedX((bestLength-1)/2+baseIndex),newFs,frameShift);%std(fftfilt(ww,x));
    f0CandidatesMap(ii,:) = crossingStructure.rawF0;
    amplitudeMap(ii,:) = crossingStructure.amplitudes;
    smoothedAMMap(ii,:) = crossingStructure.smoothedAM;
    smoothedFMMap(ii,:) = crossingStructure.smoothedFM;
    smoothedSMMap(ii,:) = crossingStructure.smoothedSM;
end;
analysisConditions.newFs = newFs;
%toc
%%
f0Structure.f0CandidatesMap = f0CandidatesMap;
f0Structure.amplitudeMap = amplitudeMap;
f0Structure.smoothedAMMap = smoothedAMMap;
f0Structure.smoothedFMMap = smoothedFMMap;
f0Structure.smoothedSMMap = smoothedSMMap;
f0Structure.temporalPositions = crossingStructure.temporalPositions;
f0Structure.fcList = fcList;
f0Structure.elapsedTime = toc;

f0Structure.analysisConditions = analysisConditions;

return;

