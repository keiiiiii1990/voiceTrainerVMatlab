function outSignalStr = ...
    generateFMspeedUpSignalRevIP(carrierFreq,initialModulationFreq,numberOfHarmonics,speedIncreaseRate,depthInSemiTone,initialPhase,signalLength,fs)
tt = 0:1/fs:signalLength;
tt = tt(:); % make it column
ge = log(speedIncreaseRate);
modulationFrequency = initialModulationFreq*exp(ge*tt);
modulationPhase = initialModulationFreq/ge*(exp(ge*tt)-1);
loMmodulationBase = sin(2*pi*modulationPhase+initialPhase)*(depthInSemiTone/12)+log2(carrierFreq);
instantFrequencyBase = 2.0.^loMmodulationBase;
instanTaneousPhase = cumsum(2*pi*instantFrequencyBase/fs);
outSignal = instanTaneousPhase*0;
for ii = 1:numberOfHarmonics
    outSignal = outSignal+sin(instanTaneousPhase*ii)/(2+ii);
end;
outSignalStr.outSignal = outSignal;
outSignalStr.instantFrequencyBase = instantFrequencyBase;
outSignalStr.modulationPhase = modulationPhase;
outSignalStr.modulationFrequency = modulationFrequency;
outSignalStr.timeBase = tt;
return;
