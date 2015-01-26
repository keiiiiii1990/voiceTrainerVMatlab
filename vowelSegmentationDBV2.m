%%  vowel segmentation based on fast F0 extractor
function segmentData = vowelSegmentationDBV2(x,fs,periodicityThreshold)
%   segmentData = vowelSegmentationDB(x,fs)

%   Originally designed by Hideki Kawahara
%   24/May/2012
%   25/Oct./2012 more reliable F0 extractor

%%
f0Structure = extractInitialF0bySymmetryICSP2013(x,fs);

%%
smoothing = 0.05; % default 50 ms
frameRate = 1000; % frame rate in Hz
w = hanning(round(smoothing*frameRate/2)*2+1);
halfLength = round(smoothing*frameRate/2);
timeAxis = f0Structure.temporalPositions;

numberOfFrames = length(f0Structure.rawPeriodicity);
smoothPeriodicity = fftfilt(w/sum(w),[f0Structure.rawPeriodicity; zeros(2*halfLength,1)]);
smoothPeriodicity = smoothPeriodicity(halfLength+(1:numberOfFrames));

vuv = double(smoothPeriodicity>periodicityThreshold); % threshold is fixed 25/Oct./2012

%%

segmentList = zeros(100,2);
segmentID = 0;
voiced = 0;
for ii = 1:numberOfFrames
    if voiced == 0 && vuv(ii) > 0.5
        voiced = 1;
        segmentID = segmentID+1;
        segmentList(segmentID,1) = timeAxis(ii);
    elseif voiced == 1 && vuv(ii) <=0.5
        voiced = 0;
        segmentList(segmentID,2) = timeAxis(ii);
    end;
end;
if voiced == 1
    segmentList(segmentID,2) = timeAxis(numberOfFrames);
end;
segmentList = segmentList(1:segmentID,:);

%%
validity = true;
numberOfSegment = segmentID;
if numberOfSegment < 5
    validity = false;
elseif numberOfSegment > 5
    [lengthList,sortedIndex] = sort(segmentList(:,2)-segmentList(:,1),'descend');
    originalList = segmentList;
    segmentID = 0;
    for ii = 1:numberOfSegment
        if ~isempty(find(sortedIndex(1:5),ii))
            segmentID = segmentID+1;
            segmentList(segmentID,:) = originalList(ii,:);
        end;
    end;
    segmentList = segmentList(1:segmentID,:);
end;
segmentData.waveform = x;
segmentData.samplingFrequency = fs;
segmentData.f0Structure = f0Structure;
segmentData.segmentList = segmentList;
segmentData.validity = validity;
%%
return;

