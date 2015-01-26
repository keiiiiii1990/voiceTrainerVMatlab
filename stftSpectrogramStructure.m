function sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms,windowShiftInms,windowType,segment,fftlSet)
%   spectral analysis for F0 test
%   sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms,windowShiftInms,windowType,segment)

%   26/Nov./2011 by Hideki Kawahara
%   10/Mar./2012 segment selection
%   08/Apr./2014 added phase spectrogram

startTime = tic;
x = x(:);
halfWindowLength = round(windowLengthInms/1000*fs/2);
fftl = 2^ceil(log2(halfWindowLength)+2);
switch nargin
    case 5
        temporalLocations = 0:windowShiftInms/1000:(length(x)-1)/fs;
    case 6
        temporalLocations = segment(1):windowShiftInms/1000:segment(2);
    case 7
        temporalLocations = segment(1):windowShiftInms/1000:segment(2);
        fftl = fftlSet;
end;
temporalLocationsInSamples = round(temporalLocations*fs);
w = feval(windowType,2*halfWindowLength+1);
baseIndex = -halfWindowLength:halfWindowLength;
nFrames = length(temporalLocations);
rawSpectrogram = zeros(fftl/2+1,nFrames);
phaseSpectrogram = rawSpectrogram;
normalizedSpectrogram = zeros(fftl/2+1,nFrames);
sgramStr.frequencyAxis = (0:fftl/2)/fftl*fs;

for ii = 1:nFrames
    tmpComplex = fft(w.*x(max(1,min(length(x),baseIndex+temporalLocationsInSamples(ii)))),fftl);
    tmpPwr = abs(tmpComplex).^2;
    rawSpectrogram(:,ii) = tmpPwr(1:fftl/2+1);
    phaseSpectrogram(:,ii) = tmpComplex(1:fftl/2+1).*exp(1i*halfWindowLength/fs*sgramStr.frequencyAxis'*2*pi);
    normalizedSpectrogram(:,ii) = rawSpectrogram(:,ii)/sum(rawSpectrogram(:,ii));
end;
maxPower = max(max(rawSpectrogram));
sg = 10*log10(rawSpectrogram/maxPower + 0.0000000000001);
sgramStr.dBspectrogram = sg;
sgramStr.rawSpectrogram = rawSpectrogram;
sgramStr.phaseSpectrogram = angle(phaseSpectrogram);
sgramStr.normalizedSpectrogram = normalizedSpectrogram;
sgramStr.temporalPositions = temporalLocations;
sgramStr.samplingFrequency = fs;
sgramStr.waveform = x;
sgramStr.windowLengthInms = windowLengthInms;
sgramStr.windowShiftInms = windowShiftInms;
sgramStr.windowType = windowType;
sgramStr.elapsedTime = toc(startTime);
end

