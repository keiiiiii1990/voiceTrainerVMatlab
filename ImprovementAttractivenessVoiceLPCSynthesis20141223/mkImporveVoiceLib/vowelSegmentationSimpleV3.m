function segmentStructure = vowelSegmentationSimpleV3(x,fs)
%   segmentStructure = vowelSegmentationSimpleV3(x,fs)

%   Designed and coded by Hideki Kawahara
%   28/Oct./2012 extend using spectral distance V3

%%
displayInd = 0;
wLengthInMs = 0.5;
wLengthHalf = round(fs*wLengthInMs/1000/2);
w = hanning(2*wLengthHalf+1);

x = x(:);
xExt = [x;zeros(wLengthHalf*2+1,1)];
xLpf = fftfilt(w/sum(w),xExt);
xLpf = xLpf(wLengthHalf+(1:length(x)));

smoothingLInMs = 10;
smoothingN = round(fs*smoothingLInMs/1000/2);
wsmooth = hanning(2*smoothingN+1);
wsmooth = wsmooth/sum(wsmooth);
wDcx = -wsmooth;
wDcx(smoothingN+1) = wDcx(smoothingN+1)+1;
cleanX = fftfilt(wDcx,[xLpf;zeros(smoothingN*2+1,1)]);
xLpf = cleanX(smoothingN+(1:length(x)));

f0Structure = extractInitialF0bySymmetryICSP2013(cleanX,fs);
%temporalPositions = f0Structure.temporalPositions;

smoothingLInMs = 200;
smoothingN = round(fs*smoothingLInMs/1000/2);
wsmooth = hanning(2*smoothingN+1);
wsmooth = wsmooth/sum(wsmooth);
xPower = fftfilt(wsmooth,[xLpf.^2;zeros(smoothingN*2+1,1)]);
xPower200 = xPower(smoothingN+(1:length(x)));

smoothingLInMs = 40;
smoothingN = round(fs*smoothingLInMs/1000/2);
wsmooth = hanning(2*smoothingN+1);
wsmooth = wsmooth/sum(wsmooth);
xPower = fftfilt(wsmooth,[xLpf.^2;zeros(smoothingN*2+1,1)]);
xPower40 = xPower(smoothingN+(1:length(x)));

gg = diff(xPower200);
hh = (gg(2:end).*gg(1:end-1)<0).*(gg(2:end)<0);
ix = (1:length(hh));
tmppeakLocations = ix(hh>0.5);
gg = diff(xPower40);
hh = (gg(2:end).*gg(1:end-1)<0).*(gg(2:end)>0);
ix = (1:length(hh));
dipLocations = ix(hh>0.5);

peakLocations = tmppeakLocations*0;
locationId = 1;
peakLocations(locationId) = tmppeakLocations(1);
for ii = 2:length(tmppeakLocations)
    if (tmppeakLocations(ii)-tmppeakLocations(ii-1))/fs > 0.11
        locationId = locationId+1;
        peakLocations(locationId) = tmppeakLocations(ii);
    else
        peakLocations(locationId) = round((tmppeakLocations(ii-1)+tmppeakLocations(ii))/2);
    end;
end;
peakLocations = peakLocations(peakLocations>0);

if length(peakLocations) < 5
    segmentStructure.validity = false;
    return;
end;
%%
numberOfSegment = length(peakLocations);
peakLevels = xPower40(peakLocations);
segmentCenterList = zeros(numberOfSegment,1);
segmentLevelList = zeros(numberOfSegment,1);
medianF0List = zeros(numberOfSegment,1);
segmentID = 0;
segmentList = zeros(numberOfSegment,2);
sdList = zeros(numberOfSegment,1);
terminalId = 1;
f0 = f0Structure.f0Raw;
for ii = 1:numberOfSegment
    segmentID = segmentID+1;
    segmentCenterList(segmentID) = peakLocations(ii);
    segmentLevelList(segmentID) = peakLevels(ii);
    segmentList(segmentID,1) = peakLocations(ii)-1;
    for jj = peakLocations(ii):-1:terminalId
        segmentList(segmentID,1) = jj;
        if xPower40(jj) < peakLevels(ii)*0.8
            break
        end;
    end;
    terminalId = peakLocations(ii)+1;
    if segmentID < numberOfSegment
        forwardTerminalId = peakLocations(ii+1);
    else
        forwardTerminalId = length(x);
    end;
    for jj = peakLocations(ii):forwardTerminalId
        segmentList(segmentID,2) = jj;
        if xPower40(jj) < peakLevels(ii)*0.85
            break
        end;
    end;
    sdList(ii) = std(f0(ceil(segmentList(ii,1)/fs*1000): ...
        floor(segmentList(ii,2)/fs*1000)));
    medianF0List(ii) = median(f0(ceil(segmentList(ii,1)/fs*1000): ...
        floor(segmentList(ii,2)/fs*1000)));
end;
%%
[sdSorted,indexSorted] = sort(sdList);
indexSortedRev = indexSorted(medianF0List(indexSorted)>65);
segmentCenterListRev = zeros(5,1);
segmentLevelListRev = zeros(5,1);
medianF0ListRev = zeros(5,1);
segmentListRev = zeros(5,2);
segmentID = 0;
for ii = 1:numberOfSegment
    if ~isempty(find(indexSortedRev(1:5)==ii))
        segmentID = segmentID+1;
        segmentCenterListRev(segmentID) = segmentCenterList(ii);
        segmentLevelListRev(segmentID) = segmentLevelList(ii);
        segmentListRev(segmentID,:) = segmentList(ii,:);
        medianF0ListRev(segmentID,:) = medianF0List(ii,:);
    end;
end;
%%

t0 = tic;
meanF0 = mean(medianF0ListRev);
sgramStr = stftSpectrogramStructure(x,fs,60,10,'nuttallwin12');
cumSgram = cumsum(sgramStr.rawSpectrogram.^2);
fxTruncated1 = sgramStr.frequencyAxis((sgramStr.frequencyAxis>meanF0/4) & (sgramStr.frequencyAxis<4000+meanF0));
sgramEnvelope1 = interp1(sgramStr.frequencyAxis(:),cumSgram,fxTruncated1(:)+meanF0/2)- ...
    interp1(sgramStr.frequencyAxis(:),cumSgram,fxTruncated1(:)-meanF0/2,'linear','extrap');
cumSgram2 = cumsum(sgramEnvelope1);
fxTruncated = sgramStr.frequencyAxis((sgramStr.frequencyAxis>meanF0) & (sgramStr.frequencyAxis<4000));
sgramEnvelope = interp1(fxTruncated1(:),cumSgram2,fxTruncated(:)+meanF0/2)- ...
    interp1(fxTruncated1(:),cumSgram2,fxTruncated(:)-meanF0/2,'linear','extrap');
toc(t0)
dBenvelope = 10*log10(sgramEnvelope);
centerEnvelope = 10*log10(sgramEnvelope(:,round(segmentCenterListRev/fs*100)));
segmentListFrame = round(segmentListRev/fs*100);
for ii = 1:5
    centerEnvelope(:,ii) = 10*log10(mean(sgramEnvelope(:,segmentListFrame(ii,1):segmentListFrame(ii,2)),2));
end;
dBdist = zeros(size(dBenvelope,2),5);
for ii = 1:size(dBenvelope,2)
    dBdist(ii,:) = std(dBenvelope(:,ii*ones(5,1))-centerEnvelope);
end;
dBdistMs = interp1(sgramStr.temporalPositions(:),dBdist,f0Structure.temporalPositions(:),'linear','extrap');

%%

segmentListRevBackup = segmentListRev;
segmentListRev = segmentListRev*0;
temporalPosition = f0Structure.temporalPositions;
f0 = f0Structure.f0Raw;
terminalID = 1;
for ii = 1:5
    for jj = ceil(segmentListRevBackup(ii,1)/fs*1000)+1:-1:terminalID
        sampleID = ceil(jj/1000*fs);
        segmentListRev(ii,1) = ceil(jj/1000*fs);
        if (xPower40(sampleID) < segmentLevelListRev(ii)*0.15) || ...
                dBdistMs(jj,ii) > 8
            break
        else
            segmentListRev(ii,1) = ceil(jj/1000*fs)+1;
        end;
    end;
    terminalID = round(segmentListRevBackup(ii,2)/fs*1000)+1;
    if ii < 5
        forwardTerminalId = round(segmentListRevBackup(ii+1,1)/fs*1000);
    else
        forwardTerminalId = length(temporalPosition);
    end;
    for jj = floor(segmentListRevBackup(ii,2)/fs*1000):forwardTerminalId
        sampleID = floor(jj/1000*fs);
        segmentListRev(ii,2) = floor(jj/1000*fs);
        if sampleID >= length(xPower40) || (xPower40(sampleID) < segmentLevelListRev(ii)*0.3) || ...
                dBdistMs(jj,ii) > 8
            break
        else
            segmentListRev(ii,2) = floor(jj/1000*fs);
        end;
    end;
end;

if prod(segmentListRev(:)) == 0
    segmentListRev = segmentListRevBackup;
    disp('Segment refinement failed!');
end;

%%
if displayInd == 1
    figure;plot(f0Structure.f0Raw);grid on;
    hold on;
    stem(segmentListRevBackup(:,1)/fs*1000,1.3*median(f0Structure.f0Raw)*ones(5,1),'g')
    stem(segmentListRevBackup(:,2)/fs*1000,1.3*median(f0Structure.f0Raw)*ones(5,1),'g')
    stem(segmentListRev(:,1)/fs*1000,1.3*median(f0Structure.f0Raw)*ones(5,1),'b')
    stem(segmentListRev(:,2)/fs*1000,1.3*median(f0Structure.f0Raw)*ones(5,1),'r')
end;

%%
segmentStructure.validity = true;
segmentStructure.segmentCenterList = segmentCenterListRev;
segmentStructure.segmentLevelList = segmentLevelListRev;
segmentStructure.segmentList = segmentListRev/fs;
segmentStructure.waveform = x;
segmentStructure.samplingFrequency = fs;
segmentStructure.medianF0List = medianF0ListRev;
segmentStructure.f0Structure = f0Structure;
segmentStructure.dBenvelope = dBenvelope;
segmentStructure.fxTruncated = fxTruncated;

return;