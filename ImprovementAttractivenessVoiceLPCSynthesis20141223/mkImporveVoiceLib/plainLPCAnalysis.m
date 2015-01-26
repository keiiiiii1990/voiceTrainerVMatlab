function lpcStructure = plainLPCAnalysis(x,fs,frameShiftInMs)
startTic = tic;

r1 = extractInitialF0bySymmetryICSP2013(x,fs);
f0 = r1.f0Raw*0+median(r1.f0Raw(r1.rawPeriodicity>0/7))*0.7;
periodicityLevel = f0*0+0.9;
%xClean = removeLF(x,fsForAnalysis,f0,periodicityLevel);
xClean = removeLF(x,fs,f0,periodicityLevel);
%r1 = extractInitialF0bySymmetryICSP2013(xClean,fsForAnalysis);
r1 = extractInitialF0bySymmetryICSP2013(xClean,fs);
rClean = cleaningF0RawTrajectoryV2(r1);
%refR1 = refineExGxOutputRevSqx(xClean,fsForAnalysis,rClean);
refR1 = refineExGxOutputRevSqx(xClean,fs,rClean);
rClean.f0 = refR1.f0Refined;
rClean.f0(refR1.rawPeriodicity<0.3) = median(rClean.f0(refR1.rawPeriodicity>0.7));
%f1 = exSpectrumTSTRAIGHTGB(xSegment,fsForAnalysis,rClean);
lpcStructure.rClean = rClean;
lpcStructure.f1byplainLPCAna = exSpectrumTSTRAIGHTGB(x,fs,rClean);

windowLengthInMs = 30;
%windowLengthInMs = 5;
preEmphasis = 0.925;
%temporalPosition = (0:frameShiftInMs/1000:(length(x)-1)/fs)';
baseIndex = (-round(windowLengthInMs/2/1000*fs):round(windowLengthInMs/2/1000*fs))';

frameShiftInMs = rClean.analysisConditions.frameShift*1000;
temporalPosition = (0:frameShiftInMs/1000:(length(x)-1)/fs)';

fftl = 2^ceil(log2(length(baseIndex))+1);
rawFAxis = (0:fftl-1)/fftl*8000;
w = blackman(length(baseIndex));
w = w/sqrt(sum(w.^2));
nOrder = 11;
signalLength = length(x);
nFrames = length(temporalPosition);
lpcMatrix = zeros(nFrames,nOrder+1);
logAreaMatrix = zeros(nFrames,nOrder+1);
errList = zeros(nFrames,1);
cepstrumList = zeros(nFrames,2);
x = [0;x(2:end)-preEmphasis*x(1:end-1)];
powerList = zeros(nFrames,1);



for ii = 1:nFrames
    xSegment = x(max(1,min(signalLength,round(temporalPosition(ii)*fs)+baseIndex)));
    %ac = real(ifft(abs(fft(xSegment.*w,fftl)).^2));
    pw = abs(fft(xSegment.*w,fftl)).^2;
    powerList(ii) = sum(pw);
    lpw = log(pw+sum(pw)/fftl/1000);
    lpwn = lpw-mean(lpw);
    theta = (0:fftl-1)/fftl*2*pi;
    %c1 = 2*sum(cos(theta(:)).*lpwn)/fftl;
    c1 = 1*sum(cos(theta(:)).*lpwn)/fftl;
    c2 = sum(cos(2*theta(:)).*lpwn)/fftl; % This may not be necessary.
    
    pwc = real(exp(lpwn-c1*cos(theta(:))-c2*cos(2*theta(:))));
    ac = real(ifft(pwc));
    [alp,err,k] = levinson(ac,nOrder);
    lpcMatrix(ii,:) = alp(:)';
    cepstrumList(ii,:) = [c1 c2];
    errList(ii) = err;
    s = ref2area(-k);
    logAreaMatrix(ii,:) = log(s(:)');
end
lpcStructure.nframe = nFrames;
lpcStructure.temporalPosition = temporalPosition;
lpcStructure.powerList = powerList;
lpcStructure.lpcMatrix = lpcMatrix;
lpcStructure.logAreaMatrix = logAreaMatrix;
lpcStructure.cepstrumList = cepstrumList;
lpcStructure.errorList = errList;
lpcStructure.windowLengthInMs = windowLengthInMs;
lpcStructure.nOrder = nOrder;
lpcStructure.rawFAxis = rawFAxis;
lpcStructure.frameShiftInMs = frameShiftInMs;
lpcStructure.fftl = fftl;
lpcStructure.preEmphasis = preEmphasis;
lpcStructure.elapsedTime = toc(startTic);
end