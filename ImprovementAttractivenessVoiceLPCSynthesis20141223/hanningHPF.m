function xClean = hanningHPF(x,fs,cutOff)
t0InSample = round(fs/cutOff);
w = hanning(2*t0InSample);
w = -w/sum(w);
w(t0InSample+1) = w(t0InSample+1)+1;
xClean = fftfilt(w,[x(:,1);zeros(length(w),1)]);
xClean = xClean(t0InSample+(1:length(x)));
return;