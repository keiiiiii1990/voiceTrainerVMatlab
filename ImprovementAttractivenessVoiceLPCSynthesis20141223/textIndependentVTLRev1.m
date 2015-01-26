function vtlStructure = textIndependentVTLRev1(x,fs,analysisConditions)
%   distanceStr = canonicalVowelTemplateBasedVTLFast2(x,fs,templateForEval,evalCondition)
%   textIndependent VTL estimation

%%   Designed and coded by Hideki Kawahara
%   30/June/2013 created
%   05/July/2013 speed up
%   07/July/2013 age and gender information
%   13/July/2013 speed up version
%   09/Nov./2013 output is in cm
%   18/Nov./2013 generalized for analysis conditions

global vtlEstimatorGlobal

if nargin < 1
    initializeVTLestimator;
    vtlStructure.windowLengthInms = vtlEstimatorGlobal.windowLengthInms;
    vtlStructure.windowShiftInms = vtlEstimatorGlobal.windowShiftInms;
    vtlStructure.windowType = vtlEstimatorGlobal.windowType;
    vtlStructure.periodicityThreshold = vtlEstimatorGlobal.periodicityThreshold;
    vtlStructure.distanceBias = vtlEstimatorGlobal.distanceBias;
    vtlStructure.magnifier = vtlEstimatorGlobal.magnifier;
    return;
elseif nargin == 2
    if isempty(vtlEstimatorGlobal)
        initializeVTLestimator;
    end;
    windowLengthInms = vtlEstimatorGlobal.windowLengthInms;
    windowShiftInms = vtlEstimatorGlobal.windowShiftInms;
    windowType = vtlEstimatorGlobal.windowType;
    periodicityThreshold = vtlEstimatorGlobal.periodicityThreshold;
    distanceBias = vtlEstimatorGlobal.distanceBias;
    magnifier = vtlEstimatorGlobal.magnifier;
elseif nargin == 3
    windowLengthInms = analysisConditions.windowLengthInms;
    windowShiftInms = analysisConditions.windowShiftInms;
    windowType = analysisConditions.windowType;
    periodicityThreshold = analysisConditions.periodicityThreshold;
    distanceBias = analysisConditions.distanceBias;
    magnifier = analysisConditions.magnifier;
else
    disp('irrelevant argument. usages are');
    disp('textIndependentVTLRev1');
    disp('vtlStructure = textIndependentVTLRev1(x,fs)');
    disp('vtlStructure = textIndependentVTLRev1(x,fs,analysisConditions)');
    vtlStructure = [];
end;
disp('Warning: This routine (textIndependentVTLRev1) uses a global variable (vtlEstimatorGlobal)');

startTic = tic;
templateForEval = vtlEstimatorGlobal.templateNForEval;
sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms,windowShiftInms,windowType);
f0Structure = extractInitialF0byHigherSymmetry(x,fs);
rawSpectrogram = sgramStr.rawSpectrogram;
f0Excerpted = interp1(f0Structure.temporalPositions,f0Structure.f0Raw, ...
    sgramStr.temporalPositions,'nearest','extrap');
periodExcerpted = interp1(f0Structure.temporalPositions, ...
    f0Structure.rawPeriodicity,sgramStr.temporalPositions,'nearest','extrap');
selectedSgram = rawSpectrogram(:,periodExcerpted>periodicityThreshold);%0.75, 0.7, 0.65
selectedF0 = f0Excerpted(periodExcerpted>periodicityThreshold);%0.75, 0.7, 0.65
%%
straightSgram = selectedSgram*0;
for ii = 1:size(selectedSgram,2);
    straightSgram(:,ii) = doubleF0Eliminator(selectedSgram(:,ii),selectedF0(ii),fs);
end;
fL = templateForEval.evaluationSpec.fL;
fH = templateForEval.evaluationSpec.fH;
fN = templateForEval.evaluationSpec.fN;
fW = templateForEval.evaluationSpec.fW;
chInOctave = templateForEval.evaluationSpec.chInOctave;
normalizedSgramStr = spectrumNormalizerFast(straightSgram,[fL fH],[fN fW],chInOctave,fs);
normalizedSgram = normalizedSgramStr.sgramNormal;
%%
vtlList = 0.55:0.14:1.6;
scoreList = zeros(length(vtlList),1);
for jj = 1:length(vtlList)
    scoreList(jj) = -integratedDistanceDDFix(sgramStr,normalizedSgram,templateForEval, ...
        vtlList(jj),distanceBias,magnifier);
end;
[~,srtIndex] = min(scoreList);
vtlMin = vtlList(max(1,srtIndex-1));
vtlMax = vtlList(min(length(scoreList),srtIndex+1));

vtlBest = fminbnd(@(vtl) -integratedDistanceDDFix(sgramStr,normalizedSgram,templateForEval, ...
    vtl,distanceBias,magnifier),vtlMin,vtlMax);
%%
vtlStructure.f0Structure = f0Structure;
vtlStructure.scoreList = scoreList;
vtlStructure.vtlList = vtlList;
vtlStructure.vtlBest = vtlBest;
vtlStructure.vtlBestInCm = vtlBest*14.66-0.998;
vtlStructure.windowLengthInms = windowLengthInms;
vtlStructure.windowShiftInms = windowShiftInms;
vtlStructure.windowType = windowType;
vtlStructure.periodicityThreshold = periodicityThreshold;
vtlStructure.distanceBias = distanceBias;
vtlStructure.magnifier = magnifier;
vtlStructure.elapsedTime = toc(startTic);
end

function initializeVTLestimator
global vtlEstimatorGlobal
gg = mfilename('fullpath');
templateName = [gg(1:strfind(gg,'textIndependentVTLRev1')-1), 'templateSet'];
tmp = load(templateName);
vtlEstimatorGlobal.templateNForEval = tmp.templateNForEval;
vtlEstimatorGlobal.windowLengthInms = 100;
vtlEstimatorGlobal.windowShiftInms = 50;
vtlEstimatorGlobal.windowType = 'blackman';
vtlEstimatorGlobal.periodicityThreshold = 0.65;
vtlEstimatorGlobal.distanceBias = 1.6;
vtlEstimatorGlobal.magnifier = 0.5;
end