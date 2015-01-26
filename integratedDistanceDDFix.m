function totalScore = integratedDistanceDDFix(sgramStr,normalizedSgram,templateForEval,vtl,bias,magn)
%   This is a test using diverse density-like integration
%   10/Nov./2013 fixed version

nFrames = size(normalizedSgram,2);
nOfSpeakers = size(templateForEval.templateAdB,2);
mapA = zeros(nFrames,nOfSpeakers);mapI = mapA;mapU = mapA;mapE = mapA;mapO = mapA;
sampledNSgramdB = 10*log10(interp1(sgramStr.frequencyAxis,normalizedSgram, ...
    templateForEval.frequencySampler/vtl));
for ii = 1:nFrames
    plobeMx = sampledNSgramdB(:,ii)*ones(1,nOfSpeakers);
    mapA(ii,:) = std(templateForEval.templateAdB-plobeMx);
    mapI(ii,:) = std(templateForEval.templateIdB-plobeMx);
    mapU(ii,:) = std(templateForEval.templateUdB-plobeMx);
    mapE(ii,:) = std(templateForEval.templateEdB-plobeMx);
    mapO(ii,:) = std(templateForEval.templateOdB-plobeMx);
end;
%bias = 2;
%magn = 2;
[nFrame,nTalker] = size(mapA);
gaussA = exp(-magn*(mapA-bias).^2);
gaussI = exp(-magn*(mapI-bias).^2);
gaussU = exp(-magn*(mapU-bias).^2);
gaussE = exp(-magn*(mapE-bias).^2);
gaussO = exp(-magn*(mapO-bias).^2);
maskA = 1-gaussA;
maskI = 1-gaussI;
maskU = 1-gaussU;
maskE = 1-gaussE;
maskO = 1-gaussO;
ddA = gaussA.*maskI.*maskU.*maskE.*maskO;
ddI = gaussI.*maskA.*maskU.*maskE.*maskO;
ddU = gaussU.*maskA.*maskI.*maskE.*maskO;
ddE = gaussE.*maskA.*maskU.*maskI.*maskO;
ddO = gaussO.*maskA.*maskU.*maskE.*maskI;
totalScore = sum(ddA(:))+sum(ddI(:))+sum(ddU(:))+sum(ddE(:))+sum(ddO(:));
totalScore = totalScore/length(ddA(:))/5;
end