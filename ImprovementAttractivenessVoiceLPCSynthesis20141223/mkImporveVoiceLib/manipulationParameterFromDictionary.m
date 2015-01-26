function manipulationStructure = manipulationParameterFromDictionary(Attractiveness,UnAttractiveness)

%example
%attractive = load('testAttExtPram.mat')
%unattractive = load('testUnAttExtPram.mat')
manipulationStructure = struct;
%AttParameterStr = Attractiveness.extParam;
%UnAttParameterStr = UnAttractiveness.extParam;
AttParameterStr = Attractiveness.exPram;
UnAttParameterStr = UnAttractiveness.exPram;

%% calculate f0 manipulation

%%%median
attMedian = AttParameterStr.f0.medianF0;
unAttMedian = UnAttParameterStr.f0.medianF0;

%medianManipulator = unAttMedian - attMedian;
medianManipulator = attMedian - unAttMedian;
medianRatioManipulator = attMedian/unAttMedian;

%%%std
stdDiff = UnAttParameterStr.f0.stdF0 - AttParameterStr.f0.stdF0;
%%%rangeF0
rangeF0Att = AttParameterStr.f0.maxF0 - AttParameterStr.f0.minF0;
rangeF0UnAtt = UnAttParameterStr.f0.maxF0 - AttParameterStr.f0.minF0;

manipulationStructure.manif0.medianManipulator = medianManipulator;
manipulationStructure.manif0.medianRatioManipulator = medianRatioManipulator;
manipulationStructure.manif0.stdDiff = stdDiff;
manipulationStructure.manif0.rangeF0Att = rangeF0Att;
manipulationStructure.manif0.rangeF0UnAtt = rangeF0UnAtt;

%% calculate rvtl manipulation

attBestRvtl = AttParameterStr.rvtl.distStr.vtlBest;
unAttBestRvtl = UnAttParameterStr.rvtl.distStr.vtlBest;

rvtlManipulator = attBestRvtl/unAttBestRvtl;
manipulationStructure.rvtl.rvtlManipulator = rvtlManipulator;
%% calculate vtaf manipulation

%%% for Japanese vowel sound DifflogVTAF 
attVtaf = AttParameterStr.vtaf.VtafStruct.lpcStructure.logAreaMatrix;
attVtafLen = length(attVtaf);
unAttVtaf = UnAttParameterStr.vtaf.VtafStruct.lpcStructure.logAreaMatrix;
unAttVtafLen = length(unAttVtaf);

if attVtafLen > unAttVtafLen
    [m,n] = size(unAttVtaf);
    newAttVtaf = zeros(m,n);
    newAttVtaf = attVtaf(m,n);
    
elseif attVtafLen < unAttVtafLen
    [m,n] = size(unAttVtaf);
    newAttVtaf = zeros(m,n);
    for ii = 1:n
        newAttVtaf(:,ii) = interp1(1:length(attVtaf),attVtaf(:,ii),1:m);
        for jj = 1:length(newAttVtaf(:,ii))
            if isnan(newAttVtaf(jj,ii))
                newAttVtaf(jj,ii) = newAttVtaf(jj-1,ii);
            end
        end
    end
    
else
    newAttVtaf = attVtaf;
end

%manipulatorLogVtaf = unAttVtaf - newAttVtaf;
manipulatorLogVtaf = newAttVtaf - unAttVtaf;
manipulationStructure.vtaf.manipulatorLogVtaf = manipulatorLogVtaf;

%%%diff median area log vtaf
diffmedianLogVTAF = UnAttParameterStr.vtaf.areaMedian - AttParameterStr.vtaf.areaMedian;
medianRatioLogVTAF = AttParameterStr.vtaf.areaMedian/UnAttParameterStr.vtaf.areaMedian;

manipulationStructure.vtaf.diffmedianLogVTAF = diffmedianLogVTAF;
manipulationStructure.vtaf.medianRatioLogVTAF = medianRatioLogVTAF;

return;