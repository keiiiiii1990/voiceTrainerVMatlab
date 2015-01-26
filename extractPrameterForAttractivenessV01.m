function exPram = extractPrameterForAttractivenessV01(x,fs)
%[x,fs] = wavread('./voice/testAttSVT01.wav');
% first : this version is for Vowel Speech
exPram = struct;

%% F0 relation parameter

f0Struct = exF0candidatesTSTRAIGHTGB(x,fs);
%medianF0 = median(f0Struct.F0(f0Struct.latentSDcoeff<1.07));
voicedF0 = f0Struct.f0(f0Struct.periodicityLevel>0.5);
unVoicedF0 = f0Struct.f0(f0Struct.periodicityLevel<0.5);

medianF0 = median(f0Struct.f0(f0Struct.periodicityLevel>0.5));
stdF0 = std(f0Struct.f0(f0Struct.periodicityLevel>0.5));
maxF0 = max(f0Struct.f0(f0Struct.periodicityLevel>0.5));
minF0 = min(f0Struct.f0(f0Struct.periodicityLevel>0.5));

exPram.f0.f0Struct = f0Struct;
exPram.f0.medianF0 = medianF0;
exPram.f0.stdF0 = stdF0;
exPram.f0.maxF0 = maxF0;
exPram.f0.minF0 = minF0;

%% vtl relation Prameter

load('templateForEval31.mat');

evalConditionCn.clipLevel = 1.6; 
evalConditionCn.pedestalLevel = -5;

y = hanningHPF(hanningHPF(hanningHPF(x,fs,70),fs,70),fs,70);
distStr = canonicalVowelTemplateBasedVTLFast2(y,fs,templateForEval31,evalConditionCn);

exPram.rvtl.evalConditionCn = evalConditionCn;
exPram.rvtl.distStr = distStr;


%% vtaf relation Parameter
ax = x;samplingFrequency = fs;
VtafStruct = struct;
VtafStruct.ax = ax;
VtafStruct.samplingFrequency = samplingFrequency;

frameLengthInMs = 25;
frameShiftInMs = 5;

responseLengthInMs = 20;
Fs = 8000;

convertedStr = convertTo8kHzSampledSignal(ax,samplingFrequency);
lpcStructure = plainLPCAnalysis(convertedStr.signal,convertedStr.samplingFrequency,frameShiftInMs);
VtafStruct.lpcStructure = lpcStructure;

f0Struct = higherSymKalmanWithTIFupdate(ax,samplingFrequency);
%f0Struct = exF0candidatesTSTRAIGHTGB(x,fs);

signalTime = (0:1/Fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);

for ii = 1:Fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/Fs*2*pi)).*(f0InSignalTime*ii<Fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/Fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*Fs):round(lpcStructure.frameShiftInMs/1000*Fs))';
responseLength = round(responseLengthInMs/1000*Fs);
fftLSynth = 2^ceil(log2(responseLength+length(baseIndex)+1));
olaIndex = (1:fftLSynth)';

for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStructure.preEmphasis*cos(theta(:))) ...
        +lpcStructure.cepstrumList(ii,1)*cos(theta(:)) ...
        +lpcStructure.cepstrumList(ii,2)*cos(2*theta(:));
    LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStructure.powerList(ii)/sqrt(f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(fftl/2:fftl) = 0;
    cepstrum(2:fftl/2) = cepstrum(2:fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = excitation(max(1,min(length(excitation),baseIndex+round(lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:responseLength);
    response = real(ifft(fft(filterResponse,fftLSynth).*fft(segment,fftLSynth)));
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*Fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*Fs)))+response;
end;

waveform = outputBuffer;
Fs = 8000;
VtafStruct.waveform = waveform;
VtafStruct.Fs = Fs;

exPram.vtaf.VtafStruct = VtafStruct;

%temporary extract

areaLen = length(exPram.vtaf.VtafStruct.lpcStructure.logAreaMatrix(1,:));
areaMedian = zeros(1,areaLen);
areaStd = zeros(1,areaLen);

areaMax = zeros(1,areaLen);
areaMaxIndex = zeros(1,areaLen);
areaMin = zeros(1,areaLen);
areaMinIndex = zeros(1,areaLen);

for ii = 1:areaLen
    areaMedian(ii) = median(exPram.vtaf.VtafStruct.lpcStructure.logAreaMatrix(:,ii));
    areaStd(ii) = std(exPram.vtaf.VtafStruct.lpcStructure.logAreaMatrix(:,ii));
    [m,ind] = max(exPram.vtaf.VtafStruct.lpcStructure.logAreaMatrix(:,ii));
    areaMax(ii) = m;
    areaMaxIndex(ii) = ind;
    [m,ind] = min(exPram.vtaf.VtafStruct.lpcStructure.logAreaMatrix(:,ii));
    areaMin(ii) = m;
    areaMinIndex(ii) = ind;
end;

exPram.vtaf.areaMedian = areaMedian;
exPram.vtaf.areaStd = areaStd;
exPram.vtaf.areaMax = areaMax;
exPram.vtaf.areaMaxIndex = areaMaxIndex;
exPram.vtaf.areaMin = areaMin;
exPram.vtaf.areaMinIndex = areaMinIndex;

return;