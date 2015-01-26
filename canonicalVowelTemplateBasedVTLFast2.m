function distanceStr = canonicalVowelTemplateBasedVTLFast2(x,fs,templateForEval,evalCondition)
%   distanceStr = canonicalVowelTemplateBasedVTLFast2(x,fs,templateForEval,evalCondition)
%   VTL estimation without labeling

%%   Designed and coded by Hideki Kawahara
%   30/June/2013 created
%   05/July/2013 speed up
%   07/July/2013 age and gender information
%   13/July/2013 speed up version

startTic = tic;
windowLengthInms = 100;
windowShiftInms = 20;
windowType = 'blackman';
sgramStr = stftSpectrogramStructure(x,fs,windowLengthInms,windowShiftInms,windowType);
f0Structure = extractInitialF0byHigherSymmetry(x,fs);
rawSpectrogram = sgramStr.rawSpectrogram;
f0Excerpted = interp1(f0Structure.temporalPositions,f0Structure.f0Raw, ...
    sgramStr.temporalPositions,'nearest','extrap');
periodExcerpted = interp1(f0Structure.temporalPositions, ...
    f0Structure.rawPeriodicity,sgramStr.temporalPositions,'nearest','extrap');
selectedSgram = rawSpectrogram(:,periodExcerpted>0.75);
selectedF0 = f0Excerpted(periodExcerpted>0.75);
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
%sampledNSgramdB = 10*log10(normalizedSgram(evalFlistIndex,:));
nFrames = size(normalizedSgram,2);
nOfSpeakers = size(templateForEval.templateAdB,2);
mapA = zeros(nFrames,nOfSpeakers);mapI = mapA;mapU = mapA;mapE = mapA;mapO = mapA;
clipLevel = evalCondition.clipLevel;
pedestalLevel = evalCondition.pedestalLevel;
vtlList = 0.55:0.14:1.6;
scoreList = zeros(length(vtlList),1);
for jj = 1:length(vtlList)
    scoreList(jj) = integratedDistance(sgramStr,normalizedSgram,templateForEval, ...
        vtlList(jj),clipLevel,pedestalLevel);%mean(logMap2(:)/5);
end;
[sortedScore,srtIndex] = min(scoreList);
vtlMin = vtlList(max(1,srtIndex-1));
vtlMax = vtlList(min(length(scoreList),srtIndex+1));

vtlBest = fminbnd(@(vtl) integratedDistance(sgramStr,normalizedSgram,templateForEval, ...
        vtl,clipLevel,pedestalLevel),vtlMin,vtlMax);
%%
distanceStr.f0Structure = f0Structure;
%distanceStr.evalCondition = evalCondition;
distanceStr.scoreList = scoreList;
distanceStr.vtlList = vtlList;
distanceStr.vtlBest = vtlBest;
distanceStr.elapsedTime = toc(startTic);
return;
