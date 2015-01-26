%%  Make vowel template with median F0
function vowelTemplate = vowelTemplateGeneratorForDBV3(segmentData)
%   vowelTemplate = vowelTemplateGeneratorForDB(segmentData)

%   Originally designed by Hideki Kawahara
%   24/May/2012
%   27/Oct./2012 extended segment with power and f0 : V2
%   28/Oct./2012 extended segment with power, f0 and spDstnc: V3

x = segmentData.waveform;
fs = segmentData.samplingFrequency;
%%
f0Structure = segmentData.f0Structure;% extractInitialF0bySymmetryICSP2013(x,fs);
%temporalPositions = f0Structure.temporalPositions;
segmentList = segmentData.segmentList;
windowLengthInms = 60;
windowShiftInms = 10;
windowType = 'nuttallwin12';

sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms, ...
    windowShiftInms,windowType,segmentList(1,:));
%%
frequencyAxis = sgramStr.frequencyAxis;
rawTemplate = zeros(length(frequencyAxis),5);
segmentMedianF0 = zeros(5,1);
segmentFrameCount = zeros(5,1);

for ii = 1:5
    sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms, ...
        windowShiftInms,windowType,segmentList(ii,:));
    segmentFrameCount(ii) = length(sgramStr.temporalPositions);
    rawTemplate(:,ii) = mean(sgramStr.rawSpectrogram,2);
    %baseFrame = (temporalPositions > segmentList(ii,1)) & ...
    %    (temporalPositions < segmentList(ii,2));
    segmentMedianF0(ii) = segmentData.medianF0List(ii);%median(f0Structure.f0Raw(baseFrame));
end;
%%
vowelTemplate.rawTemplate = rawTemplate;
vowelTemplate.segmentMedianF0 = segmentMedianF0;
vowelTemplate.segmentFrameCount = segmentFrameCount;
vowelTemplate.frequencyAxis = frequencyAxis;
vowelTemplate.f0Structure = f0Structure;
vowelTemplate.segmentData = segmentData;
return;

