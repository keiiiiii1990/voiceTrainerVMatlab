function convertedStr = convertTo8kHzSampledSignal(x,fs)
%   convert signal to 8kHz sampled one (quick and dirty implementation)
%   by Hideki Kawahara
tx = (0:1/8000:(length(x)-1)/fs)';
l2p5ms = round(0.0025*fs);
baseIndex = -l2p5ms:l2p5ms;
fftl = 4096;
fx = (0:fftl-1)'/fftl*fs;
fxx = fx;
fxx(fx>fs/2) = fxx(fx>fs/2)-fs;
filtDesign = abs(fxx)<4000;
baseFIR = fftshift(real(ifft(filtDesign)));
trimmedFIR = baseFIR(baseIndex+fftl/2+1).*hanning(length(baseIndex));
filterdSignal = fftfilt(trimmedFIR,x);
convertedStr.signal = interp1(((0:length(x)-1)-l2p5ms)/fs,filterdSignal,tx,'linear',0);
convertedStr.temporalPosition = tx;
convertedStr.FIRfilter = trimmedFIR;
convertedStr.originalSamplingFrequency = fs;
convertedStr.samplingFrequency = 8000;
convertedStr.filterdSignal = filterdSignal;
end