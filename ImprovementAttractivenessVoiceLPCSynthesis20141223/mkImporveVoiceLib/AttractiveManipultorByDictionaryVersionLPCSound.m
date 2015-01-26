function manipulateSoundStructure = AttractiveManipultorByDictionaryVersionLPCSound(x,fs,dictionary)

%[x,fs] = wavread('./voice/m1705.wav');
%dictionary = load('hmAttractiveness.mat')
%load(dictionary)
mStr = dictionary.manipulationStructure;
manipulateSoundStructure = struct;

%% vtaf relation Parameter
ax = x;samplingFrequency = fs;
lpcStruct = struct;
lpcStruct.ax = ax;
lpcStruct.samplingFrequency = samplingFrequency;

frameLengthInMs = 25;
frameShiftInMs = 5;

responseLengthInMs = 20;
Fs = 8000;

lpcStruct.convertedStr = convertTo8kHzSampledSignal(ax,samplingFrequency);
lpcStructure = plainLPCAnalysis(lpcStruct.convertedStr.signal,lpcStruct.convertedStr.samplingFrequency,frameShiftInMs);
lpcStruct.lpcStructure = lpcStructure;

lpcStruct.f0Struct = higherSymKalmanWithTIFupdate(ax,samplingFrequency);
%f0Struct = exF0candidatesTSTRAIGHTGB(x,fs);

lpcStruct.signalTime = (0:1/Fs:lpcStruct.f0Struct.temporalPositions(end))';
lpcStruct.excitation = lpcStruct.signalTime*0;
lpcStruct.f0InSignalTime = interp1(lpcStruct.f0Struct.temporalPositions,lpcStruct.f0Struct.F0,lpcStruct.signalTime,'linear',50);

for ii = 1:Fs/40/2
    lpcStruct.excitation = lpcStruct.excitation+cos(cumsum(ii*lpcStruct.f0InSignalTime/Fs*2*pi)).*(lpcStruct.f0InSignalTime*ii<Fs/2);
end;

lpcStruct.outputBuffer = lpcStruct.excitation*0;
lpcStruct.theta = lpcStruct.lpcStructure.rawFAxis/Fs*2*pi;
lpcStruct.fftl = length(lpcStruct.theta);
lpcStruct.baseIndex = (-round(lpcStruct.lpcStructure.frameShiftInMs/1000*Fs):round(lpcStruct.lpcStructure.frameShiftInMs/1000*Fs))';
lpcStruct.responseLength = round(responseLengthInMs/1000*Fs);
lpcStruct.fftLSynth = 2^ceil(log2(lpcStruct.responseLength+length(lpcStruct.baseIndex)+1));
lpcStruct.olaIndex = (1:lpcStruct.fftLSynth)';

for ii = 1:length(lpcStruct.lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStruct.lpcStructure.preEmphasis*cos(lpcStruct.theta(:))) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,1)*cos(lpcStruct.theta(:)) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,2)*cos(2*lpcStruct.theta(:));
    LPCspectrumLn = -2*log(abs(fft(lpcStruct.lpcStructure.lpcMatrix(ii,:)',lpcStruct.fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStruct.lpcStructure.powerList(ii)/sqrt(lpcStruct.f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(lpcStruct.fftl/2:lpcStruct.fftl) = 0;
    cepstrum(2:lpcStruct.fftl/2) = cepstrum(2:lpcStruct.fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = lpcStruct.excitation(max(1,min(length(lpcStruct.excitation),lpcStruct.baseIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:lpcStruct.responseLength);
    response = real(ifft(fft(filterResponse,lpcStruct.fftLSynth).*fft(segment,lpcStruct.fftLSynth)));
    lpcStruct.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs))) = ...
        lpcStruct.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs)))+response;
end;

waveform = lpcStruct.outputBuffer;
Fs = 8000;
lpcStruct.waveform = waveform;
lpcStruct.Fs = Fs;

manipulateSoundStructure.lpcStruct = lpcStruct;
%% f0 manipulaiton
%%{
mf0 = mStr.manif0.medianManipulator;
manipulateSoundStructure.f0.mf0 = mf0;
mRatef0 = mStr.manif0.medianRatioManipulator;
manipulateSoundStructure.f0.mRatef0 = mRatef0;

manipulateSoundStructure.f0.f0InSignalTime = interp1(lpcStruct.f0Struct.temporalPositions,lpcStruct.f0Struct.F0,lpcStruct.signalTime,'linear',50)*mRatef0;
%manipulateSoundStructure.f0.f0InSignalTime = interp1(lpcStruct.f0Struct.temporalPositions,lpcStruct.f0Struct.F0,lpcStruct.signalTime,'linear',50);
%manipulateSoundStructure.f0.f0InSignalTime = manipulateSoundStructure.f0.f0InSignalTime - mf0;

manipulateSoundStructure.f0.excitation = lpcStruct.excitation;
for ii = 1:fs/40/2
    manipulateSoundStructure.f0.excitation = manipulateSoundStructure.f0.excitation + cos(cumsum(ii*manipulateSoundStructure.f0.f0InSignalTime/Fs*2*pi)).*(manipulateSoundStructure.f0.f0InSignalTime*ii<Fs/2);
end;

manipulateSoundStructure.f0.outputBuffer = manipulateSoundStructure.f0.excitation*0;

for ii = 1:length(lpcStruct.lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStruct.lpcStructure.preEmphasis*cos(lpcStruct.theta(:))) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,1)*cos(lpcStruct.theta(:)) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,2)*cos(2*lpcStruct.theta(:));
    LPCspectrumLn = -2*log(abs(fft(lpcStruct.lpcStructure.lpcMatrix(ii,:)',lpcStruct.fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStruct.lpcStructure.powerList(ii)/sqrt(lpcStruct.f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(lpcStruct.fftl/2:lpcStruct.fftl) = 0;
    cepstrum(2:lpcStruct.fftl/2) = cepstrum(2:lpcStruct.fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = manipulateSoundStructure.f0.excitation(max(1,min(length(manipulateSoundStructure.f0.excitation),lpcStruct.baseIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:lpcStruct.responseLength);
    response = real(ifft(fft(filterResponse,lpcStruct.fftLSynth).*fft(segment,lpcStruct.fftLSynth)));
    manipulateSoundStructure.f0.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs))) = ...
        manipulateSoundStructure.f0.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs)))+response;
end;
%}

%% rvtl manipulation
%%{
mrvtl = mStr.rvtl.rvtlManipulator;
manipulateSoundStructure.vtl.mrvtl = mrvtl;
manipulateSoundStructure.vtl.frequencyAxis = (0:lpcStruct.fftl-1)'/lpcStruct.fftl*8000;

manipulateSoundStructure.vtl.outputBuffer = lpcStruct.excitation*0;

for ii = 1:length(lpcStruct.lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStruct.lpcStructure.preEmphasis*cos(lpcStruct.theta(:))) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,1)*cos(lpcStruct.theta(:)) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,2)*cos(2*lpcStruct.theta(:));
    LPCspectrumLn = -2*log(abs(fft(lpcStruct.lpcStructure.lpcMatrix(ii,:)',lpcStruct.fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedHalfLPCspectrum = interp1(manipulateSoundStructure.vtl.frequencyAxis,fixedLPCspectrum, ...
        mrvtl*manipulateSoundStructure.vtl.frequencyAxis(1:lpcStruct.fftl/2));
    fixedLPCspectrum = [fixedHalfLPCspectrum;fixedHalfLPCspectrum(end-1:-1:2)];
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStruct.lpcStructure.powerList(ii)/sqrt(lpcStruct.f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(lpcStruct.fftl/2:lpcStruct.fftl) = 0;
    cepstrum(2:lpcStruct.fftl/2) = cepstrum(2:lpcStruct.fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = lpcStruct.excitation(max(1,min(length(lpcStruct.excitation),lpcStruct.baseIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:lpcStruct.responseLength);
    response = real(ifft(fft(filterResponse,lpcStruct.fftLSynth).*fft(segment,lpcStruct.fftLSynth)));
    manipulateSoundStructure.vtl.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs))) = ...
        manipulateSoundStructure.vtl.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs)))+response;
end;
%}
%% vtaf manipulation

%%{
mLogVtaf = mStr.vtaf.manipulatorLogVtaf;
nSection = length(mLogVtaf(1,:));

if length(mLogVtaf(:,1)) > length(lpcStruct.lpcStructure.lpcMatrix(:,1))
    newMLogVtaf = mLogVtaf(1:length(lpcStruct.lpcStructure.lpcMatrix(:,1)),:);

elseif length(mLogVtaf(:,1)) < length(lpcStruct.lpcStructure.lpcMatrix(:,1))
    newMLogVtaf = zeros(length(lpcStruct.lpcStructure.lpcMatrix(:,1)),nSection);
    for ii = 1:nSection
        newMLogVtaf(:,ii) = interp1(1:length(mLogVtaf(:,ii)),mLogVtaf(:,ii),1:length(lpcStruct.lpcStructure.lpcMatrix(:,1)));
    end
else
    newMLogVtaf = mLogVtaf;
end

manipulateSoundStructure.vtaf.mLogVtaf = newMLogVtaf;
manipulateSoundStructure.vtaf.outputBuffer = lpcStruct.excitation*0;

for ii = 1:length(lpcStruct.lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStruct.lpcStructure.preEmphasis*cos(lpcStruct.theta(:))) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,1)*cos(lpcStruct.theta(:)) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,2)*cos(2*lpcStruct.theta(:));
    logArea = lpcStruct.lpcStructure.logAreaMatrix(ii,:);
    newLogArea = logArea+newMLogVtaf(ii,:)/4;
    newRef = area2ref(exp(newLogArea));                                     %%newRef is reflection coefficient
    newalp = k2alp(newRef);                                                 %%newalp is predictor
    LPCspectrumLn = -2*log(abs(fft(newalp(:),lpcStruct.fftl)));                       %%calculate LPC spectrum
    %LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStruct.lpcStructure.powerList(ii)/sqrt(lpcStruct.f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(lpcStruct.fftl/2:lpcStruct.fftl) = 0;
    cepstrum(2:lpcStruct.fftl/2) = cepstrum(2:lpcStruct.fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = lpcStruct.excitation(max(1,min(length(lpcStruct.excitation),lpcStruct.baseIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:lpcStruct.responseLength);
    response = real(ifft(fft(filterResponse,lpcStruct.fftLSynth).*fft(segment,lpcStruct.fftLSynth)));
    manipulateSoundStructure.vtaf.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs))) = ...
        manipulateSoundStructure.vtaf.outputBuffer(min(length(lpcStruct.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs)))+response;
end;
%}
%% all parameter manipultion

manipulateSoundStructure.all.outputBuffer = lpcStruct.excitation*0;
for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStruct.lpcStructure.preEmphasis*cos(lpcStruct.theta(:))) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,1)*cos(lpcStruct.theta(:)) ...
        +lpcStruct.lpcStructure.cepstrumList(ii,2)*cos(2*lpcStruct.theta(:));
    logArea = lpcStruct.lpcStructure.logAreaMatrix(ii,:);
    newLogArea = logArea+manipulateSoundStructure.vtaf.mLogVtaf(ii,:)/4;
    newRef = area2ref(exp(newLogArea));                                     %%newRef is reflection coefficient
    newalp = k2alp(newRef);                                                 %%newalp is predictor
    LPCspectrumLn = -2*log(abs(fft(newalp(:),lpcStruct.fftl)));
    %LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedHalfLPCspectrum = interp1(manipulateSoundStructure.vtl.frequencyAxis,fixedLPCspectrum, ...
        manipulateSoundStructure.vtl.mrvtl*manipulateSoundStructure.vtl.frequencyAxis(1:lpcStruct.fftl/2));
    fixedLPCspectrum = [fixedHalfLPCspectrum;fixedHalfLPCspectrum(end-1:-1:2)];
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStruct.lpcStructure.powerList(ii)/sqrt(lpcStruct.f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(lpcStruct.fftl/2:lpcStruct.fftl) = 0;
    cepstrum(2:lpcStruct.fftl/2) = cepstrum(2:lpcStruct.fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = manipulateSoundStructure.f0.excitation(max(1,min(length(manipulateSoundStructure.f0.excitation),lpcStruct.baseIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:lpcStruct.responseLength);
    response = real(ifft(fft(filterResponse,lpcStruct.fftLSynth).*fft(segment,lpcStruct.fftLSynth)));
    manipulateSoundStructure.all.outputBuffer(min(length(manipulateSoundStructure.f0.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs))) = ...
        manipulateSoundStructure.all.outputBuffer(min(length(manipulateSoundStructure.f0.excitation),lpcStruct.olaIndex+round(lpcStruct.lpcStructure.temporalPosition(ii)*Fs)))+response;
end;

return;