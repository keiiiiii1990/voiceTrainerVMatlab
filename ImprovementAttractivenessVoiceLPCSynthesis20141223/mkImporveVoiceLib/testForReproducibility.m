%%  Test script for reproducibility of distance analysis
%   by Hideki Kawahara
%   30/Oct./2012
%   01/Nov./2012 explicit seed assignment
%   03/Nov./2012 reploducibility rest

%   typical setting
dataDirectory = 'vowels20121028T185523OK';
% testType = 'flozen'; % 'random' or 'full'
testType = 'full'; % 'random' or 'full'
maxTalker = 100;
baseWidth = 2500;
smoothingWidth = 300;
boundary.upper = 5000;
boundary.lower = 600;
displayOn = 0;
seedForRand = 5;

%%

distanceStructureBase = pairedDistanceAnalysis(dataDirectory, ...
    'full',maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);

%%

maxTalkerList = [25 50 100];
iterations = 10;
testType = 'random';

testResults = struct;
standardErrorList = zeros(iterations,length(maxTalkerList));
sDofEstimatedVTLR = zeros(iterations,length(maxTalkerList));

testID = 0;
for ii = 1:length(maxTalkerList)
    maxTalker = maxTalkerList(ii);
    for jj = 1:iterations
        distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
            testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);
        
        standardErrorList(jj,ii) = distanceStructure.standardError;
        sDofEstimatedVTLR(jj,ii) = std(distanceStructure.estimatedVTLR(:));
        testID = testID+1;
        testResults(testID).distanceStructure = distanceStructure;
        testResults(testID).iterationID = jj;
        testResults(testID).maxTalker = maxTalker;
    end;
end;

%%

figure;plot(sort(standardErrorList(:,1)),(1:10)/10);grid on;
hold all
plot(sort(standardErrorList(:,2)),(1:10)/10);grid on;
plot(sort(standardErrorList(:,3)),(1:10)/10);grid on;

figure;plot(sort(sDofEstimatedVTLR(:,1)),(1:10)/10);grid on;
hold all
plot(sort(sDofEstimatedVTLR(:,2)),(1:10)/10);grid on;
plot(sort(sDofEstimatedVTLR(:,3)),(1:10)/10);grid on;

figure;plot(sort(standardErrorList(:,1)./sDofEstimatedVTLR(:,1)),(1:10)/10);grid on;
hold all
plot(sort(standardErrorList(:,2)./sDofEstimatedVTLR(:,2)),(1:10)/10);grid on;
plot(sort(standardErrorList(:,3)./sDofEstimatedVTLR(:,3)),(1:10)/10);grid on;

%%

figure;
errorlocations = -0.2:0.01:0.2;
for ii = 1:testID
    vtlError = testResults(ii).distanceStructure.VTLRatioMatrix(:) - ...
        testResults(ii).distanceStructure.estimatedVTLR(:);
    probabilityList = (1:length(vtlError))/length(vtlError);
    plot(sort(vtlError),probabilityList);grid on;
    hold all;
end;
hold off;

%%

figure;
plot(maxTalkerList,mean(standardErrorList./sDofEstimatedVTLR),'o-', ...
    'linewidth',3);grid on;
hold on;
plot(maxTalkerList,mean(standardErrorList./sDofEstimatedVTLR)+ ...
    std(standardErrorList./sDofEstimatedVTLR)/2,'o-');
plot(maxTalkerList,mean(standardErrorList./sDofEstimatedVTLR)- ...
    std(standardErrorList./sDofEstimatedVTLR)/2,'o-');
plot(maxTalkerList,min(standardErrorList./sDofEstimatedVTLR),'*');
plot(maxTalkerList,max(standardErrorList./sDofEstimatedVTLR),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('relative error');
outName = ['repRelError' num2str(now,30) '.eps'];
print('-deps',outName);

%%

maxTalkerListX = [maxTalkerList 385];
relError385 = distanceStructureBase.standardError/distanceStructureBase.SDofVTLR;
figure;
semilogx(maxTalkerListX,[mean(standardErrorList./sDofEstimatedVTLR) relError385],'o-', ...
    'linewidth',3);grid on;
hold on;
semilogx(maxTalkerList,mean(standardErrorList./sDofEstimatedVTLR)+ ...
    std(standardErrorList./sDofEstimatedVTLR)/2,'o-');
semilogx(maxTalkerList,mean(standardErrorList./sDofEstimatedVTLR)- ...
    std(standardErrorList./sDofEstimatedVTLR)/2,'o-');
semilogx(maxTalkerList,min(standardErrorList./sDofEstimatedVTLR),'*');
semilogx(maxTalkerList,max(standardErrorList./sDofEstimatedVTLR),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('relative error');
outName = ['repRelErrorL' num2str(now,30) '.eps'];
print('-deps',outName);

%%

figure;
plot(maxTalkerList,mean(standardErrorList),'o-', ...
    'linewidth',3);grid on;
hold on;
plot(maxTalkerList,mean(standardErrorList)+ ...
    std(standardErrorList)/2,'o-');
plot(maxTalkerList,mean(standardErrorList)- ...
    std(standardErrorList)/2,'o-');
plot(maxTalkerList,min(standardErrorList),'*');
plot(maxTalkerList,max(standardErrorList),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('rms error of VTL ratio');
outName = ['repVTLRlError' num2str(now,30) '.eps'];
print('-deps',outName);

%%

stdError385 = distanceStructureBase.standardError;
figure;
semilogx(maxTalkerListX,[mean(standardErrorList)  stdError385],'o-', ...
    'linewidth',3);grid on;
hold on;
semilogx(maxTalkerList,mean(standardErrorList)+ ...
    std(standardErrorList)/2,'o-');
semilogx(maxTalkerList,mean(standardErrorList)- ...
    std(standardErrorList)/2,'o-');
semilogx(maxTalkerList,min(standardErrorList),'*');
semilogx(maxTalkerList,max(standardErrorList),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('rms error of VTL ratio');
outName = ['repVTLRlErrorL' num2str(now,30) '.eps'];
print('-deps',outName);

%%

figure;
plot(maxTalkerList,mean(sDofEstimatedVTLR),'o-', ...
    'linewidth',3);grid on;
hold on;
plot(maxTalkerList,mean(sDofEstimatedVTLR)+ ...
    std(sDofEstimatedVTLR)/2,'o-');
plot(maxTalkerList,mean(sDofEstimatedVTLR)- ...
    std(sDofEstimatedVTLR)/2,'o-');
plot(maxTalkerList,min(sDofEstimatedVTLR),'*');
plot(maxTalkerList,max(sDofEstimatedVTLR),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('SD of VTL ratio');
outName = ['repVTLrSD' num2str(now,30) '.eps'];
print('-deps',outName);

%%

SD385 = distanceStructureBase.SDofVTLR;
figure;
semilogx(maxTalkerListX,[mean(sDofEstimatedVTLR) SD385],'o-', ...
    'linewidth',3);grid on;
hold on;
semilogx(maxTalkerList,mean(sDofEstimatedVTLR)+ ...
    std(sDofEstimatedVTLR)/2,'o-');
semilogx(maxTalkerList,mean(sDofEstimatedVTLR)- ...
    std(sDofEstimatedVTLR)/2,'o-');
semilogx(maxTalkerList,min(sDofEstimatedVTLR),'*');
semilogx(maxTalkerList,max(sDofEstimatedVTLR),'*');
set(gca,'fontsize',15);
legend('mean','mean+sd/2','mean-sd/2','minimum','maximum');
xlabel('number of talkers');
ylabel('SD of VTL ratio');
outName = ['repVTLrSDL' num2str(now,30) '.eps'];
print('-deps',outName);

%%

figure;
plot(sDofEstimatedVTLR(:,3),standardErrorList(:,3),'o','linewidth',2);
hold on;
plot(0.125+0.013*[-1 1],0.02+0.001*[-1 1],'--');
axis([0.11 0.14 0.019 0.021]);
grid on;
set(gca,'fontsize',16);
xlabel('SD of VTL ratio');
ylabel('SD of VTL ratio error')
print -deps VTLsdVSestimationError.eps

%%
figure;
elapsedTimeList = [250,158,64,17];
loglog([maxTalkerList 385],elapsedTimeList,'o-');
hold on;
loglog([maxTalkerList 385],[maxTalkerList 385]./elapsedTimeList*60,'o-', ...
    'linewidth',3);
grid on;
set(gca,'fontsize',15);
xlabel('number of talkers');
ylabel('count or time (s)');
legend('count per minutes','processing time');
outName = ['speedTest' num2str(now,30) '.eps'];
print('-deps',outName);



