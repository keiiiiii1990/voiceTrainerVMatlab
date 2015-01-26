function vtlStructure = textIndependentVTLKalman(x,fs,analysisConditions)
%   vtlStructure = textIndependentVTLKalman(x,fs,analysisConditions)
%   textIndependent VTL estimation

%%   Designed and coded by Hideki Kawahara
%   30/June/2013 created
%   05/July/2013 speed up
%   07/July/2013 age and gender information
%   13/July/2013 speed up version
%   09/Nov./2013 output is in cm
%   18/Nov./2013 generalized for analysis conditions

global vtlEstimatorGlobal

if isempty(vtlEstimatorGlobal)
    disp('Global bariable: vtlEstimatorGlobal is empty. Initialization is invoked.v');
    initializeVTLestimator;
end;

if nargin < 1
    initializeVTLestimator;
    vtlStructure.windowLengthInms = vtlEstimatorGlobal.windowLengthInms;
    vtlStructure.windowShiftInms = vtlEstimatorGlobal.windowShiftInms;
    vtlStructure.windowType = vtlEstimatorGlobal.windowType;
    vtlStructure.periodicityThreshold = vtlEstimatorGlobal.periodicityThreshold;
    vtlStructure.distanceBias = vtlEstimatorGlobal.distanceBias;
    vtlStructure.magnifier = vtlEstimatorGlobal.magnifier;
    disp('Warning: This routine (textIndependentVTLRev1) uses a global variable (vtlEstimatorGlobal)');
    return;
elseif nargin == 2
    if isempty(vtlEstimatorGlobal)
        initializeVTLestimator;
        disp('Warning: This routine (textIndependentVTLRev1) uses a global variable (vtlEstimatorGlobal)');
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

startTic = tic;
templateForEval = vtlEstimatorGlobal.templateNForEval;
sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms,windowShiftInms,windowType);
f0Structure = higherSymKalmanInitial(x,fs);
rawSpectrogram = sgramStr.rawSpectrogram;
reliableIndex = f0Structure.latentSDcoeff < periodicityThreshold;
f0Excerpted = interp1(f0Structure.temporalPositions(reliableIndex),f0Structure.latentF0(reliableIndex), ...
    sgramStr.temporalPositions,'nearest','extrap');
latentSDExcerpted = interp1(f0Structure.temporalPositions(reliableIndex), ...
    f0Structure.latentSDcoeff(reliableIndex),sgramStr.temporalPositions,'nearest','extrap');
selectedSgram = rawSpectrogram(:,latentSDExcerpted<periodicityThreshold);
selectedF0 = f0Excerpted(latentSDExcerpted<periodicityThreshold);
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
%vtlList = 0.5:0.1:1.4;
vtlList = 0.65*2.0.^((0:8)/8*log2(1.4/0.65));
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
vtlStructure.evaluationSpec = templateForEval.evaluationSpec;
vtlStructure.elapsedTime = toc(startTic);
end

function initializeVTLestimator
global vtlEstimatorGlobal
gg = mfilename('fullpath');
if ~isempty(strfind(gg,'textIndependentVTLKalman'))
    templateName = [gg(1:strfind(gg,'textIndependentVTLKalman')-1), 'templateSet'];
    tmp = load(templateName);
else
    disp('Vowel template file: templateSet.mat');
    clear global vtlEstimatorGlobal
    return;
end;
vtlEstimatorGlobal.templateNForEval = tmp.templateNForEval;
vtlEstimatorGlobal.windowLengthInms = 100;
vtlEstimatorGlobal.windowShiftInms = 50;
vtlEstimatorGlobal.windowType = 'blackman';
vtlEstimatorGlobal.periodicityThreshold = 1.04;
vtlEstimatorGlobal.distanceBias = 1.6;
vtlEstimatorGlobal.magnifier = 0.5;
end