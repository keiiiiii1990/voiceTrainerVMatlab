function mixedMeasure = individualToMixedMeasure(r,normalizeConstant,weight,alphaP,betaP)

AMn = r.smoothedAMMap/normalizeConstant(1);
FMn = r.smoothedFMMap/normalizeConstant(2);
SMn = r.smoothedSMMap/normalizeConstant(3);
mixedMeasure = exp(-alphaP*(weight(1)*AMn.^betaP+weight(2)*FMn.^betaP+weight(3)*SMn.^betaP).^(1/betaP));

return;