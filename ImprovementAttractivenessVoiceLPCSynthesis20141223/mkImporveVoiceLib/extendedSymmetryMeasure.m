%%  Test for extended zero crossing engine
function intervalStructure = extendedSymmetryMeasure(xs,fs,frameShift)
%   intervalStructure = extendedSymmetryMeasure(xs,fs,frameShift)
%   xs  : input signal under study
%   fs  : sampling frequency (Hz)
%   frameShift  : frame shift period (s)

%   designed and coded by Hideki Kawahara
%   17/August/2011
%   09/Oct./2011
%   10/Oct./2011
%   12/Oct./2011 debug
%   16/Oct./2011 optimization based on profiling
%   07/Apr./2012 revised for Interspeech 2012
%   14/June/2012 revised for level information
%   27/Aug./2012 higher order asymmetry (mid point)
%   01/Sep./2012 re-designed
%   02/Sep./2012 small revision
%   23/Sep./2012 minor bug fix
%   12/Sep./2013 bug fix
%%

xs = xs(:)-mean(xs(:));
diff1 = [diff(xs); 0];
diff2 = [0; diff(xs)];
outLocations = (0:frameShift:(length(xs)-1)/fs)'; % 1 ms frame 

temporalBin = (1:length(xs))';
extremes = (diff1.*diff2)<0; 
ipoint = 0;
tmpextremeIndices = temporalBin(extremes);
extremeIndices = tmpextremeIndices*0;
lastIndex = -10000;
for jj = 1:length(tmpextremeIndices)
    if lastIndex ~= tmpextremeIndices(jj)
        ipoint = ipoint+1;
        extremeIndices(ipoint) = tmpextremeIndices(jj);
        lastIndex = tmpextremeIndices(jj);
    else
        disp(['duplicate at:' num2str(tmpextremeIndices(jj))])
    end;
end;
extremeIndices = extremeIndices(1:ipoint);

if length(extremeIndices) > 7  % minor bug after local smoothing 23/Sept/2012
    extremeIndices = extremeIndices(2:end-1);
    extremeLocations = ...
        diff2(extremeIndices)./(diff2(extremeIndices)-diff2(extremeIndices+1))-0.5+extremeIndices;
    if sum(isnan(extremeLocations))>0;extremeLocations = extremeIndices;disp('NaN is detected');end;
    extremeValues = ...
        -((xs(extremeIndices+1)-xs(extremeIndices-1)).^2)./(8*(xs(extremeIndices-1)+xs(extremeIndices+1)-2*xs(extremeIndices))) ...
        + xs(extremeIndices);
    extremeIntervals = diff(extremeLocations)/fs; 
    extremeAmplitudes = abs(diff(extremeValues));
    extremeT0Intervals = extremeIntervals(1:end-1)+extremeIntervals(2:end);
    extremeMeanAmplitudes = extremeAmplitudes(1:end-1)+extremeAmplitudes(2:end);
    midLocations = (extremeLocations(1:end-1)+extremeLocations(2:end))/2;
    midLevels = (extremeValues(1:end-1)+extremeValues(2:end))/2;
    midValues = interp1q(temporalBin,xs,midLocations);
    relativeSM = abs(midLevels-midValues)./extremeAmplitudes;
    relativeFM = abs(diff(extremeIntervals))./extremeT0Intervals;
    relativeAM = abs(diff(extremeAmplitudes))./extremeMeanAmplitudes;
    segmentF0 = 1.0./extremeT0Intervals;
    extremeTemporalLocations = extremeLocations(2:end-1)/fs;
    midTemporalLocation = midLocations/fs;
    smoothedAM = (relativeAM(1:end-2)+relativeAM(2:end-1)+relativeAM(3:end))/3;
    smoothedFM = (relativeFM(1:end-2)+relativeFM(2:end-1)+relativeFM(3:end))/3;
    smoothedSM = (relativeSM(1:end-2)+relativeSM(2:end-1)+relativeSM(3:end))/3;
    extremeTemporalLocations2 = extremeLocations(3:end-2)/fs;
    midTemporalLocation2 = midLocations(2:end-1)/fs;
    
    rawF0 = interp1q([0;extremeTemporalLocations(:);length(xs)/fs],[segmentF0(1);segmentF0(:);segmentF0(end)],outLocations);
    relativeFMinterpolated = ...
        interp1q([0;extremeTemporalLocations(:);length(xs)/fs],[relativeFM(1);relativeFM(:);relativeFM(end)],outLocations);
    relativeAMinterpolated = ...
        interp1q([0;extremeTemporalLocations(:);length(xs)/fs],[relativeAM(1);relativeAM(:);relativeAM(end)],outLocations);
    relativeSMinterpolated = ...
        interp1q([0;midTemporalLocation(:);length(xs)/fs],[relativeSM(1);relativeSM(:);relativeSM(end)],outLocations);
    amplitudes = ...
        interp1q([0;extremeTemporalLocations(:);length(xs)/fs],[extremeMeanAmplitudes(1);extremeMeanAmplitudes(:);extremeMeanAmplitudes(end)],outLocations);

    smoothedFMinterpolated = ...
        interp1q([0;extremeTemporalLocations2(:);length(xs)/fs],[smoothedFM(1);smoothedFM(:);smoothedFM(end)],outLocations);
    smoothedAMinterpolated = ...
        interp1q([0;extremeTemporalLocations2(:);length(xs)/fs],[smoothedAM(1);smoothedAM(:);smoothedAM(end)],outLocations);
    smoothedSMinterpolated = ...
        interp1q([0;midTemporalLocation2(:);length(xs)/fs],[smoothedSM(1);smoothedSM(:);smoothedSM(end)],outLocations);
    

    intervalStructure.rawF0 = rawF0;
    intervalStructure.relativeFM = relativeFMinterpolated;
    intervalStructure.relativeAM = relativeAMinterpolated;
    intervalStructure.relativeSM = relativeSMinterpolated;
    intervalStructure.smoothedFM = smoothedFMinterpolated;
    intervalStructure.smoothedAM = smoothedAMinterpolated;
    intervalStructure.smoothedSM = smoothedSMinterpolated;
    intervalStructure.amplitudes = amplitudes;
    intervalStructure.temporalPositions = outLocations;
else
    intervalStructure.rawF0 = outLocations*0+1;
    intervalStructure.relativeFM = outLocations*0+1;
    intervalStructure.relativeAM = outLocations*0+1;
    intervalStructure.relativeSM = outLocations*0+1;
    intervalStructure.smoothedFM = outLocations*0+1;
    intervalStructure.smoothedAM = outLocations*0+1;
    intervalStructure.smoothedSM = outLocations*0+1;
    intervalStructure.amplitudes = outLocations*0+0.0001;
    intervalStructure.temporalPositions = outLocations;
end;
%%
intervalStructure.temporalPositions = outLocations;
return;

