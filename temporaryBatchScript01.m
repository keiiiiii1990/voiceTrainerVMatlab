%%
%[ax,samplingFrequency] = wavread('./voice/Attractiveness_Vowel02_All_SYtr.wav');
%[ax,samplingFrequency] = wavread('./voice/Attractiveness_Sentence_01_SY.wav');
%[ax,samplingFrequency] = wavread('./voice/Attractiveness_Vowel01_all_SYtr.wav');
%[ax,samplingFrequency] = wavread('./voice/UnAttractiveness_Vowel01_all_SYtr.wav');
[ax,samplingFrequency] = wavread('./voice/testAttSVT01.wav');
%[ax,samplingFrequency] = wavread('./voice/testAttSVT01.wav');

AttractiveStruct = struct;
AttractiveStruct.ax = ax;
AttractiveStruct.samplingFrequency = samplingFrequency;
%[uax,samplingFrequency] = wavread('');

frameLengthInMs = 25;
frameShiftInMs = 5;

responseLengthInMs = 20;
fs = 8000;

convertedStr = convertTo8kHzSampledSignal(ax,samplingFrequency);
lpcStructure = plainLPCAnalysis(convertedStr.signal,convertedStr.samplingFrequency,frameShiftInMs);
AttractiveStruct.lpcStructure = lpcStructure;

%lpcStructure = myGUIdata.lpcStructure;
f0Struct = higherSymKalmanWithTIFupdate(ax,samplingFrequency);
%sgramStr = myGUIdata.sgramStr;
signalTime = (0:1/fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);

for ii = 1:fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/fs*2*pi)).*(f0InSignalTime*ii<fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*fs):round(lpcStructure.frameShiftInMs/1000*fs))';
responseLength = round(responseLengthInMs/1000*fs);
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
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs)))+response;
end;

x = outputBuffer;
fs = 8000;

AttractiveStruct.waveform = x;
AttractiveStruct.fs = fs;


%logarea mean
[m,n] = size(AttractiveStruct.lpcStructure.logAreaMatrix);
logAreaAtt = zeros(m+1,n);
logAreaAtt(1:459,:) = AttractiveStruct.lpcStructure.logAreaMatrix;
logAreaAtt(end,:) = AttractiveStruct.lpcStructure.logAreaMatrix(end,:);
logAreaAtt(end-3:end,:)
logAreaAttVowel = zeros(5,12);
cnt = 0;
for ii = 1:5
    for jj = 1:length(logAreaAtt(:,1))/5
        cnt = cnt + 1;
        logAreaAttVowel(ii,:) = logAreaAttVowel(ii,:) + logAreaAtt(cnt,:);
    end
    logAreaAttVowel(ii,:) = logAreaAttVowel(ii,:)/92;
end



%% unattractive vtaf

[ax,samplingFrequency] = wavread('./voice/testUnAttSVT01.wav');

UnAttractiveStruct = struct;
UnAttractiveStruct.ax = ax;
UnAttractiveStruct.samplingFrequency = samplingFrequency;
%[uax,samplingFrequency] = wavread('');

frameLengthInMs = 25;
frameShiftInMs = 5;

responseLengthInMs = 20;
fs = 8000;

convertedStr = convertTo8kHzSampledSignal(ax,samplingFrequency);
lpcStructure = plainLPCAnalysis(convertedStr.signal,convertedStr.samplingFrequency,frameShiftInMs);
UnAttractiveStruct.lpcStructure = lpcStructure;

%lpcStructure = myGUIdata.lpcStructure;
f0Struct = higherSymKalmanWithTIFupdate(ax,samplingFrequency);
%sgramStr = myGUIdata.sgramStr;
signalTime = (0:1/fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);

for ii = 1:fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/fs*2*pi)).*(f0InSignalTime*ii<fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*fs):round(lpcStructure.frameShiftInMs/1000*fs))';
responseLength = round(responseLengthInMs/1000*fs);
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
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs)))+response;
end;

x = outputBuffer;
fs = 8000;

UnAttractiveStruct.waveform = x;
UnAttractiveStruct.fs = fs;

%
[m,n] = size(UnAttractiveStruct.lpcStructure.logAreaMatrix);
logAreaUnAtt = zeros(m+3,n);
logAreaUnAtt(1:m,:) = UnAttractiveStruct.lpcStructure.logAreaMatrix;
logAreaUnAtt(end,:) = UnAttractiveStruct.lpcStructure.logAreaMatrix(end,:);
logAreaUnAttVowel = zeros(5,12);
cnt = 0;
for ii = 1:5
    for jj = 1:length(logAreaUnAtt(:,1))/5
        cnt = cnt + 1;
        logAreaUnAttVowel(ii,:) = logAreaUnAttVowel(ii,:) + logAreaUnAtt(cnt,:);
    end
    logAreaUnAttVowel(ii,:) = logAreaUnAttVowel(ii,:)/(length(logAreaUnAtt(:,1))/5);
end


%% vtl estimation Ver 20140713
%reference:VTLbyCanonicalDB.pdf
%reference:/Users/amlab/MATLABFORSTRAIGHT/mfiles3/progs20130713

%[x,fs] = wavread('./voice/testAttSVT01.wav');
%[x,fs] = wavread('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141118/voice/Attractiveness_Vowel01_all_SYtr.wav');
[x,fs] = wavread('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141118/voice/UnAttractiveness_Vowel01_all_SYtr.wav');
%Attractiveness_Vowel01_all_SYtr.wav
%on referece folder
load('templateForEval31.mat');

evalConditionCn.clipLevel = 1.6; 
evalConditionCn.pedestalLevel = -5;


y = hanningHPF(hanningHPF(hanningHPF(x,fs,70),fs,70),fs,70);
distStr = canonicalVowelTemplateBasedVTLFast2(y,fs,templateForEval31,evalConditionCn);


%% f0 relation parameter extraction

%[x,fs] = wavread('./voice/testUnAttSVT01.wav');
[x,fs] = wavread('./voice/testAttSVT01.wav');

%reference :sampleVTmanipulatorGUIR2.m
f0Struct = higherSymKalmanWithTIFupdate(x,fs);
medianF0 = median(f0Struct.F0(f0Struct.latentSDcoeff<1.07));


%% make vowel template memoScript

%/Users/amlab/MATLABFORSTRAIGHT/mfiles/で実行が必要．

%これらのファイルを/Users/amlab/MATLABFORSTRAIGHT/mfiles/へ
[x,fs] = wavread('Attractiveness_Vowel02_All_SYtr.wav');
vs = vowelSegmentationSimpleV3(x,fs);
save('testAttS01.mat','vs');

load('testAtt01.mat')
[x,fs] = wavread('Attractiveness_Vowel02_All_SYtr.wav');

segmentData = vs;
vowelTemplate = vowelTemplateGeneratorForDBV3(segmentData);
vowelTemplate.dataName = 'Attractiveness_Vowel02_All_SYtr.wav';
vowelTemplate.signal = x;
vowelTemplate.samplingFrequency = fs;
vowelTemplate.gender = 'm';
vowelTemplate.age = 23;
vowelTemplate.ID = 2301;
save('testAttVT01.mat','vowelTemplate')

%2kaime
load('testAtt01.mat')
[x,fs] = wavread('Attractiveness_Vowel02_All_SYtr.wav');

segmentData = vs;
vowelTemplate = vowelTemplateGeneratorForDBV3(segmentData);
vowelTemplate.dataName = 'Attractiveness_Vowel02_All_SYtr.wav';
vowelTemplate.signal = x;
vowelTemplate.samplingFrequency = fs;
vowelTemplate.gender = 'm';
vowelTemplate.age = 23;
vowelTemplate.ID = 2301;
save('testAttVT01.mat','vowelTemplate')






%% trucate vowel segmentation
%ref:memoOnTemplate-4.pdf
%ref:/Users/amlab/MATLABFORSTRAIGHT/mfiles/vowelSegmentationSimpleV3.m

%load vowel template MAT file

%value = load('./voice/testAttVT01.mat');
%value = load('./voice/testAttSVT01.mat');
value = load('./voice/testUnAttSVT01.mat');
%value.vowelTemplate = Tvalue.vs;
saveWavFilename = './voice/testUnAttSVT01.wav';

fs = value.vowelTemplate.segmentData.samplingFrequency;
y = value.vowelTemplate.segmentData.waveform;

vowelSpecInfo = struct;
for ii = 1:5
    % choice segment range
    yStart = value.vowelTemplate.segmentData.segmentList(ii, 1);
    yEnd = value.vowelTemplate.segmentData.segmentList(ii, 2);
    yUse = y(round(yStart*fs):round(yEnd*fs));
    yUse2 = yUse.*hamming(length(yUse));
    yUse3 = zeros(1, length(yUse2));
    yUse3(1) = 0;
    for l = 2:length(yUse3)
        yUse3(l) = yUse(l) - 0.98*yUse(l-1); % preEmphasis
    end;
    switch ii
        case 1
            vowelSpecInfo.original.waveform.a = yUse2;
        case 2
            vowelSpecInfo.original.waveform.i = yUse2;
        case 3
            vowelSpecInfo.original.waveform.u = yUse2;
        case 4
            vowelSpecInfo.original.waveform.e = yUse;
        case 5
            vowelSpecInfo.original.waveform.o = yUse2;
    end;
end;

vowelSpecInfo.original.samplingFrequency = fs;

alength = length(vowelSpecInfo.original.waveform.a);
ilength = length(vowelSpecInfo.original.waveform.i);
ulength = length(vowelSpecInfo.original.waveform.u);
elength = length(vowelSpecInfo.original.waveform.e);
olength = length(vowelSpecInfo.original.waveform.o);
trvowelWaveform = zeros(alength+ilength+ulength+elength+olength,1);
nextend = alength;
trvowelWaveform(1:alength) = vowelSpecInfo.original.waveform.a;

prenextend = nextend;
nextend = nextend + ilength;
trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.i;

prenextend = nextend;
nextend = nextend + ulength;
trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.u;

prenextend = nextend;
nextend = nextend + elength;
trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.e;

prenextend = nextend;
nextend = nextend + olength;
trvowelWaveform(prenextend+1:end) = vowelSpecInfo.original.waveform.o;
audiowrite(saveWavFilename,trvowelWaveform,vowelSpecInfo.original.samplingFrequency);


