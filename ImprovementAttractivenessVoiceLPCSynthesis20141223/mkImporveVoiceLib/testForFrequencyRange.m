%%  Test script for reproducibility of distance analysis
%   by Hideki Kawahara
%   30/Oct./2012
%   01/Nov./2012 explicit seed assignment
%   03/Nov./2012 reploducibility test
%   04/Nov./2012 frequency range test

%   typical setting
dataDirectory = 'vowels20121028T185523OK';
testType = 'flozen'; % 'random' or 'full'
maxTalker = 100;
baseWidth = 2500;
smoothingWidth = 300;
%boundary.upper = 5000;
%boundary.lower = 600;
displayOn = 0;
seedForRand = 5;

%%

upperBoundaryList = [3000 4000 5000 6000 7000];
lowerBoundaryList = [300 400 500 600 700];

testResults = struct;
testID = 0;
standardErrorList = zeros(length(upperBoundaryList),length(lowerBoundaryList));
sDofEstimatedVTLR = zeros(length(upperBoundaryList),length(lowerBoundaryList));
for ii = 1:length(lowerBoundaryList)
    boundary.lower = lowerBoundaryList(ii);
    for jj = 1:length(upperBoundaryList)
        boundary.upper = upperBoundaryList(jj);
        distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
            testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);
        standardErrorList(jj,ii) = distanceStructure.standardError;
        sDofEstimatedVTLR(jj,ii) = std(distanceStructure.estimatedVTLR(:));
        testID = testID+1;
        testResults(testID).distanceStructure = distanceStructure;
        testResults(testID).upperBoundary = boundary.upper;
        testResults(testID).lowerBoundary = boundary.lower;
        testResults(testID).maxTalker = maxTalker;
        testResults(testID).baseWidth = baseWidth;
        testResults(testID).smoothingWidth = smoothingWidth;
        testResults(testID).seedForRand = seedForRand;
    end;
end;

save boudaryTest1.mat testResults
%%
figure;
imagesc([lowerBoundaryList(1) lowerBoundaryList(end)], ...
    [upperBoundaryList(1) upperBoundaryList(end)],standardErrorList);
colorbar
axis('xy');
axis('square')
set(gca,'fontsize',16)
xlabel('lower boundary (Hz)');
ylabel('upper boundary (Hz)');
title(['VTLR error basew:' num2str(baseWidth) '(Hz) smthw:' num2str(smoothingWidth) '(Hz)']);
outName = ['VTLRErrorMapB' num2str(now,30) '.eps'];
print('-depsc',outName);

figure;
imagesc([lowerBoundaryList(1) lowerBoundaryList(end)], ...
    [upperBoundaryList(1) upperBoundaryList(end)],sDofEstimatedVTLR(1:5,1:5));
axis('xy');
colorbar
axis('square')
set(gca,'fontsize',16)
xlabel('lower boundary (Hz)');
ylabel('upper boundary (Hz)');
title(['VTLR SD basew:' num2str(baseWidth) '(Hz) smthw:' num2str(smoothingWidth) '(Hz)']);
outName = ['VTLRSDMapB' num2str(now,30) '.eps'];
print('-depsc',outName);

figure;
imagesc([lowerBoundaryList(1) lowerBoundaryList(end)], ...
    [upperBoundaryList(1) upperBoundaryList(end)],standardErrorList./sDofEstimatedVTLR);
axis('xy');
colorbar
axis('square')
set(gca,'fontsize',16)
xlabel('lower boundary (Hz)');
ylabel('upper boundary (Hz)');
title(['rel error basew:' num2str(baseWidth) '(Hz) smthw:' num2str(smoothingWidth) '(Hz)']);
outName = ['relErrorMapB' num2str(now,30) '.eps'];
print('-depsc',outName);

%%

boundary.upper = 4000;
boundary.lower = 600;

baseWidthList = [1000 1400 2000 2800 4000];
smoothingWidthList = [150 210 300 425 600];

testResults2 = struct;
testID = 0;
standardErrorList = zeros(length(smoothingWidthList),length(baseWidthList));
sDofEstimatedVTLR = zeros(length(smoothingWidthList),length(baseWidthList));
for ii = 1:length(baseWidthList)
    baseWidth = baseWidthList(ii);
    for jj = 1:length(upperBoundaryList)
        smoothingWidth = smoothingWidthList(jj);
        distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
            testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);
        standardErrorList(jj,ii) = distanceStructure.standardError;
        sDofEstimatedVTLR(jj,ii) = std(distanceStructure.estimatedVTLR(:));
        testID = testID+1;
        testResults2(testID).distanceStructure = distanceStructure;
        testResults2(testID).baseWidth = baseWidth;
        testResults2(testID).smoothingWidth = smoothingWidth;
        testResults2(testID).maxTalker = maxTalker;
        testResults2(testID).baseWidth = baseWidth;
        testResults2(testID).smoothingWidth = smoothingWidth;
        testResults2(testID).seedForRand = seedForRand;
    end;
end;

save smoothingTest1.mat testResults2

%%

figure;
imagesc([baseWidthList(1) baseWidthList(end)], ...
    [smoothingWidthList(1) smoothingWidthList(end)],standardErrorList);
colorbar
axis('xy');
axis('square')
set(gca,'fontsize',16)
xlabel('trend remover (Hz)');
ylabel('smoother (Hz)');
title(['VTLR error LB:' num2str(boundary.lower) '(Hz) UB:' num2str(boundary.upper) '(Hz)']);
outName = ['VTLRErrorMapS' num2str(now,30) '.eps'];
print('-depsc',outName);

figure;
imagesc([baseWidthList(1) baseWidthList(end)], ...
    [smoothingWidthList(1) smoothingWidthList(end)],sDofEstimatedVTLR(1:5,1:5));
axis('xy');
colorbar
axis('square')
set(gca,'fontsize',16)
xlabel('trend remover (Hz)');
ylabel('smoother (Hz)');
title(['VTLR SD LB:' num2str(boundary.lower) '(Hz) UB:' num2str(boundary.upper) '(Hz)']);
outName = ['VTLRSDMapS' num2str(now,30) '.eps'];
print('-depsc',outName);

figure;
imagesc([baseWidthList(1) baseWidthList(end)], ...
    [smoothingWidthList(1) smoothingWidthList(end)],standardErrorList./sDofEstimatedVTLR);
axis('xy');
colorbar
axis('square')
set(gca,'fontsize',16)
xlabel('trend remover (Hz)');
ylabel('smoother (Hz)');
title(['rel error LB:' num2str(boundary.lower) '(Hz) UB:' num2str(boundary.upper) '(Hz)']);
outName = ['relErrorMapS' num2str(now,30) '.eps'];
print('-depsc',outName);

%%  tentative best setting

testType = 'full'; % 'flozen', 'random' or 'full'
maxTalker = 1000;
baseWidth = 2500;
smoothingWidth = 300;
boundary.upper = 4000;
boundary.lower = 600;
displayOn = 1;
seedForRand = 5;

distanceStructureBase = pairedDistanceAnalysis(dataDirectory, ...
    'full',maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);

%%
genderCharStrct = distanceStructureBase.templateNameList;
maleList = zeros(size(genderCharStrct,2),2);
femaleList = zeros(size(genderCharStrct,2),2);
femaleID = 0;
maleID = 0;
for ii = 1:size(genderCharStrct,2)
    switch genderCharStrct(ii).name(1)
        case 'f'
            femaleID = femaleID+1;
            femaleList(femaleID,1) = ii;
            femaleList(femaleID,2) = distanceStructureBase.vtl(ii);
        case 'm'
            maleID = maleID+1;
            maleList(maleID,1) = ii;
            maleList(maleID,2) = distanceStructureBase.vtl(ii);
    end;
end;

femaleList = femaleList(1:femaleID,:);
maleList = maleList(1:maleID,:);
figure;
plot(sort(femaleList(:,2)),(1:femaleID)/femaleID,'r','linewidth',2);
hold on;
plot(sort(maleList(:,2)),(1:maleID)/maleID,'b','linewidth',2);
grid on;
set(gca,'fontsize',15);
xlabel('relative VTL');
ylabel('probability');
legend('female','male');
outName = ['VTLdist' datestr(now,30) '.eps'];
print('-depsc',outName);




