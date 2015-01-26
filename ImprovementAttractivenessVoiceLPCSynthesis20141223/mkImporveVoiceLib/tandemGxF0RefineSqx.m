function sourceStructure = tandemGxF0RefineSqx(x,fs,currentTime,f0initial,numberOfHarmonics,windowStretch,numberOfIntererence)
%   F0 refinement using generalized TANDEM-based instantaneous frequency calculation
%   sourceStructure = tandemGxF0Refine(x,fs,currentTime,f0initial,numberOfHarmonics,windowStretch,numberOfIntererence)

%   Originally designed for ICASSP2011 by Hideki Kawahara
%   23/Oct./2010 initial version
%   18/Oct./2011 revised for ExZx
%   21/Oct./2011 generalization
%   23/Oct./2011 revise
%   13/Oct./2011 weighting is square root (x for test)

%numberOfIntererence = 1;
if isnan(f0initial)
    f0initial = 100;
    %disp(currentTime);
end;
if f0initial < 32
    f0initial = 32; % safe guard
end;
if f0initial > 2000
    f0initial = 2000;
end;
%halfWindowLength = ceil(1.3*6/(numberOfIntererence+1)*fs/f0initial/2);
halfWindowLength = ceil(windowStretch*fs/f0initial/2);
numberOfWindows = numberOfIntererence*2;
centerLocationIndex = (0:2*numberOfIntererence-1)';
centerLocationList = ((2*centerLocationIndex+1)/(4*numberOfIntererence)-1/2)/f0initial;
windowLengthInTime = (2*halfWindowLength+1)/fs;
baseIndex = (-halfWindowLength:halfWindowLength)';
baseTime = baseIndex/fs;
fftl = 2^ceil(log2((halfWindowLength*2+1))+1);
fx = ((0:fftl-1)/fftl*fs)';
tandemSpectrum = zeros(fftl,1);
numeratorI = zeros(fftl,1);
for ii = 1:numberOfWindows
    preTime = currentTime+baseTime+centerLocationList(ii);
    preIndexRaw = ceil(preTime*fs);
    preIndexTime = (preIndexRaw-1)/fs;
    preWindowTime = preIndexTime-currentTime-centerLocationList(ii);
    mainWindow = 0.42+0.5*cos(2*pi*preWindowTime/windowLengthInTime)+...
        0.08*cos(4*pi*preWindowTime/windowLengthInTime);
    preDW = -(diff([0;mainWindow])+diff([mainWindow;0]))/2;
    preIndex = max(1,min(length(x),preIndexRaw));
    preSpectrum = fft(x(preIndex).*mainWindow,fftl);
    preDSpectrum = fft(x(preIndex).*preDW,fftl);
    numeratorI = numeratorI+real(preSpectrum).*imag(preDSpectrum)-imag(preSpectrum).*real(preDSpectrum);
    tandemSpectrum = tandemSpectrum+abs(preSpectrum).^2;
end;
instFreq2 = fx+numeratorI./tandemSpectrum*fs/2/pi;

trimIndex = (1:2)';
indexListTrim = round(f0initial/(fs/fftl)*trimIndex)+1;
fixpList = instFreq2(indexListTrim); 
powerList = tandemSpectrum(indexListTrim);
f0initial = sum(powerList.*(fixpList./trimIndex))/sum(powerList);

if f0initial < 32
    f0initial = 32;
end;
if f0initial > 1200
    f0initial = 1200;
end;

trimIndex = (1:numberOfHarmonics)';
indexListTrim = round(f0initial/(fs/fftl)*trimIndex)+1;
fixpList = instFreq2(indexListTrim); 
powerList = sqrt(tandemSpectrum(indexListTrim));
%meanF0 = sum(powerList.*(fixpList./trimIndex))/sum(powerList);
meanF0 = sum(powerList.*(fixpList))/sum(powerList.*trimIndex);
sourceStructure.preSpectrum = preSpectrum;
sourceStructure.instFreq2 = instFreq2;
sourceStructure.fx = fx;
sourceStructure.f0 = meanF0;








