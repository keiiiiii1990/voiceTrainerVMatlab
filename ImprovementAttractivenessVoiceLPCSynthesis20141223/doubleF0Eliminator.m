function spectrumWOinterference = doubleF0Eliminator(spectrumSlice,f0,fs)
fftl = (size(spectrumSlice,1)-1)*2;
fx = (0:fftl/2)/fftl*fs;
cumulatedSpectrum = cumsum(spectrumSlice);
tmp = (interp1(fx,cumulatedSpectrum,min(fs/2,fx+f0/2))-...
    interp1(fx,cumulatedSpectrum,max(0,fx-f0/2)))/f0;
cum2nd = cumsum(tmp);
spectrumWOinterference = (interp1(fx,cum2nd,min(fs/2,fx+f0/2))-...
    interp1(fx,cum2nd,max(0,fx-f0/2)))/f0;
return;