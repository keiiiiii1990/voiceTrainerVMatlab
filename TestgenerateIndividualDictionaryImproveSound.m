function manipulateSoundStructure = TestgenerateIndividualDictionaryImproveSound(x,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use

% reference :
% generateIndivisualCombinationSoundV2(x,fs,dictionary,savefileName,fileSaveDirName).m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[x,fs] = wavread('./voice/m1705.wav');
%dictionary = load('hmAttractiveness.mat')
%load(dictionary)
%mStr = dictionary.manipulationStructure;
manipulateSoundStructure = struct;

%fileSaveDirName = '20141203eachParameterHM';
%savefileName = 'm1705DicHM';
%mkdir(fileSaveDirName)

%% vtaf relation Parameter

%soundOrgOnlyFileName = [fileSaveDirName '/' 'Original' savefileName '.wav'];

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

%% temprary calcurate distance unattractiveness voice for make Improve voice GUI
%%%calculate vtl 
y = hanningHPF(hanningHPF(hanningHPF(x,fs,70),fs,70),fs,70);
manipulateSoundStructure.tmpPrm.currenctVtlStr = textIndependentVTLKalman(y,fs);
 
vtlDis = manipulateSoundStructure.tmpPrm.currenctVtlStr.vtlBestInCm * (0.4 - 0.1);
medianF0Dis = median(manipulateSoundStructure.lpcStruct.f0Struct.F0(manipulateSoundStructure.lpcStruct.f0Struct.latentSDcoeff>1.0)) * (1-0.4);
manipulateSoundStructure.tmpPrm.tmpDistance = vtlDis + medianF0Dis;

return;