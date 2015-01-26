function pBasedWeight = periodicityIndex2VbasedWeight2(f0PPStr)
%   simpler version
%   by Hideki Kawahara
%   11/Aug./2013
%   14/Aug./2013 revised
%   18/Aug./2013 fixed minor fragility

ticInit = tic;
%biasCoefficient = 0.65; % this is the center
biasCoefficient = 0.75; % biased more
nChannel = length(f0PPStr.fcList);
nFrame = size(f0PPStr.updatedPeriodicityMap,2);
weightBuffer = zeros(nChannel,nFrame);
normalizedWeightBuffer = weightBuffer;
weightBufferTmp = zeros(nChannel,1);
nInOctave = f0PPStr.analysisConditions.nInOctave;
fcIdBuffer = nInOctave*log2(biasCoefficient*f0PPStr.f0CandidatesMap/f0PPStr.fcList(1))+1;
variationFcId = zeros(nChannel,1);
sigmaEstMap = periodicityIdx2sigmaInCent(f0PPStr.updatedPeriodicityMap);
for kk = 1:nFrame
    weightBufferTmp = weightBufferTmp*0;
    for ii = 1:nChannel
        fcId = fcIdBuffer(ii,kk);
        if (fcId>2.6) && (fcId<nChannel-1.6)
            weightBufferTmp(round(fcId)) = weightBufferTmp(round(fcId))+1/sigmaEstMap(ii,kk);
        end;
    end;
    weightBuffer(:,kk) = weightBufferTmp;
    normalizedWeightBuffer(:,kk) = weightBufferTmp/sum(weightBufferTmp);
end;
f0cent = 1200*log2(f0PPStr.f0CandidatesMap/440);
averagedF0cent = sum(f0cent.*normalizedWeightBuffer);
averagedF0 = 440*2.0.^(averagedF0cent/1200);
averagedFcId = nInOctave*log2(biasCoefficient*averagedF0/f0PPStr.fcList(1))+1;
for kk = 1:nFrame
    variationFcId(kk) = sum(normalizedWeightBuffer(:,kk).*(fcIdBuffer(:,kk)-averagedFcId(kk)).^2);
end;
averagedF0Fcid = max(1,min(nChannel,round(nInOctave*log2(biasCoefficient*averagedF0/f0PPStr.fcList(1))+1)));
pBasedWeight.averagedF0 = averagedF0;
pBasedWeight.variationFcIdInCent2 = variationFcId*200^2;
for kk = 1:nFrame
    pBasedWeight.totalSigma(kk) = pBasedWeight.variationFcIdInCent2(kk)+sigmaEstMap(averagedF0Fcid(kk),kk);
end;
pBasedWeight.weightBuffer = weightBuffer;
pBasedWeight.normalizedWeightBuffer = normalizedWeightBuffer;
pBasedWeight.sigmaEstMap = sigmaEstMap;
pBasedWeight.elapsedTime = toc(ticInit);
return;