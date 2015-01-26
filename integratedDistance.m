function totalScore = integratedDistance(sgramStr,normalizedSgram,templateForEval9,vtl,clipLevel,pedestalLevel)
nFrames = size(normalizedSgram,2);
nOfSpeakers = size(templateForEval9.templateAdB,2);
mapA = zeros(nFrames,nOfSpeakers);mapI = mapA;mapU = mapA;mapE = mapA;mapO = mapA;
sampledNSgramdB = 10*log10(interp1(sgramStr.frequencyAxis,normalizedSgram, ...
    templateForEval9.frequencySampler/vtl));
for ii = 1:nFrames
    plobeMx = sampledNSgramdB(:,ii)*ones(1,nOfSpeakers);
    mapA(ii,:) = std(templateForEval9.templateAdB-plobeMx);
    mapI(ii,:) = std(templateForEval9.templateIdB-plobeMx);
    mapU(ii,:) = std(templateForEval9.templateUdB-plobeMx);
    mapE(ii,:) = std(templateForEval9.templateEdB-plobeMx);
    mapO(ii,:) = std(templateForEval9.templateOdB-plobeMx);
end;
logMap2 = max(pedestalLevel,min(clipLevel,log(mapA)))+ ...
    max(pedestalLevel,min(clipLevel,log(mapI)))+ ...
    max(pedestalLevel,min(clipLevel,log(mapU)))+ ...
    max(pedestalLevel,min(clipLevel,log(mapE)))+ ...
    max(pedestalLevel,min(clipLevel,log(mapO)));
totalScore =  mean(logMap2(:)/5);
return;