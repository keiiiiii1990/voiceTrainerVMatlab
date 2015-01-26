function normalizedStr = spectrumNormalizerFast(straightSgram,boundary,smoothers,chPerOct,fs)
fL = boundary(1);
fH = boundary(2);
fN = smoothers(1);
fW = smoothers(2);
fList = fL*2.0.^(0:1/chPerOct:log2(fH/fL));
nFrames = size(straightSgram,2);
fftl = (size(straightSgram,1)-1)*2;
fx = (0:fftl-1)'/fftl*fs;
fxx = fx;
fxx(fxx>fs/2) = fxx(fxx>fs/2)-fs;
smootherN = (0.5+0.5*cos(pi*fxx/fN)).*(abs(fxx/fN)<1);
smootherW = (0.5+0.5*cos(pi*fxx/fW)).*(abs(fxx/fW)<1);
lifterN = (real(fft(smootherN/sum(smootherN))));
lifterW = (real(fft(smootherW/sum(smootherW))));
doubleSgram = [straightSgram;straightSgram(end-1:-1:2,:)];
cepstrumGram = real(fft(log(doubleSgram)));
cepstrumNormal = zeros(fftl,nFrames);
for ii = 1:nFrames
    cepstrumNormal(:,ii) = (lifterN-lifterW).*cepstrumGram(:,ii);
end;
sgramNormal = exp(real(ifft(cepstrumNormal)));
normalizedStr.sgramNormal = sgramNormal(1:fftl/2+1,:);
normalizedStr.fAxis = fx(1:fftl/2+1);
normalizedStr.fList = fList;
return;