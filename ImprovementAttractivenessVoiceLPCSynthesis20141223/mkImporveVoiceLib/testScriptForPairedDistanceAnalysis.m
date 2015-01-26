%%  Sample script for distance analysis
%   by Hideki Kawahara
%   30/Oct./2012
%   01/Nov./2012 explicit seed assignment

%   typical setting
dataDirectory = 'vowels20121028T185523OK';
testType = 'flozen'; % 'random' or 'full'
maxTalker = 385;
baseWidth = 2500;
smoothingWidth = 300;
boundary.upper = 5000;
boundary.lower = 600;
displayOn = 0;
seedForRand = 5;

%%  test example for baseWidth

baseWidwhList = 1500*2.0.^(-1:1/4:1);
standardErrorList = zeros(length(baseWidwhList),1);
sDofEstimatedVTLR = zeros(length(baseWidwhList),1);

for ii = 1:length(baseWidwhList)
    baseWidth = baseWidwhList(ii);
    distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
        testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);
    standardErrorList(ii) = distanceStructure.standardError;
    sDofEstimatedVTLR(ii) = std(distanceStructure.estimatedVTLR(:));
end;

%%
figure;
plot(baseWidwhList,standardErrorList,'o-')
grid on;
set(gca,'fontsize',15)
xlabel('trend remover width (Hz)');
ylabel('standard error in VTL ratio');
title([testType ' mxtlk:' num2str(maxTalker) ' smwth:' num2str(smoothingWidth) ...
    '(Hz) Uf:' num2str(boundary.upper) '(Hz) Lf:' num2str(boundary.lower) '(Hz)']);
outName = ['vtlStdErr' num2str(now,30) '.eps'];
print('-deps',outName);

figure;
plot(baseWidwhList,standardErrorList./sDofEstimatedVTLR,'o-')
grid on;
set(gca,'fontsize',15)
xlabel('smoothing width (Hz)');
ylabel('relative standard error of VTL');
title([testType ' mxtlk:' num2str(maxTalker) ' smwth:' num2str(smoothingWidth) ...
    '(Hz) Uf:' num2str(boundary.upper) '(Hz) Lf:' num2str(boundary.lower) '(Hz)']);
outName = ['vtlStdErrBs' num2str(now,30) '.eps'];
print('-deps',outName);


%%

% baseWidth = 2000;
% smoothingWidthList = [100 200 300 400 500];
% standardErrorList = zeros(length(smoothingWidthList),1);
% sDofEstimatedVTLR = zeros(length(smoothingWidthList),1);
% 
% for ii = 1:length(smoothingWidthList)
%     smoothingWidth = smoothingWidthList(ii);
%     distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
%         testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,displayOn);
%     standardErrorList(ii) = distanceStructure.standardError;
%     sDofEstimatedVTLR(ii) = std(distanceStructure.estimatedVTLR(:));
% end;

%%
% figure;
% plot(smoothingWidthList,standardErrorList,'o-')
% grid on;
% set(gca,'fontsize',15)
% xlabel('trend remover width (Hz)');
% ylabel('standard error in VTL ratio');
% title([testType ' mxtlk:' num2str(maxTalker) ' bswth:' num2str(baseWidth) ...
%     '(Hz) Uf:' num2str(boundary.upper) '(Hz) Lf:' num2str(boundary.lower) '(Hz)']);
% outName = ['vtlStdSmRaw' num2str(now,30) '.eps'];
% print('-deps',outName);
% 
% figure;
% plot(smoothingWidthList,standardErrorList./sDofEstimatedVTLR,'o-')
% grid on;
% set(gca,'fontsize',15)
% xlabel('smoothing width (Hz)');
% ylabel('relative standard error of VTL');
% title([testType ' mxtlk:' num2str(maxTalker) ' bswth:' num2str(baseWidth) ...
%     '(Hz) Uf:' num2str(boundary.upper) '(Hz) Lf:' num2str(boundary.lower) '(Hz)']);
% outName = ['vtlStdErrSm' num2str(now,30) '.eps'];
% print('-deps',outName);
% 
% %%
% 
% baseWidth = 1700;
% smoothingWidth = 250;
% distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
%     testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,1);
% 
% %%
% maxTalker = 1000;
% 
% baseWidth = 1700;
% smoothingWidth = 250;
% testType = 'full';
% distanceStructure = pairedDistanceAnalysis(dataDirectory, ...
%     testType,maxTalker,baseWidth,smoothingWidth,boundary,seedForRand,1);
% 
% %%
% 
% ser = (distanceStructure.VTLRatioMatrix(:)-distanceStructure.estimatedVTLR(:));
% figure;plot(sort(ser),(1:length(ser))/length(ser));grid on;
% set(gca,'fontsize',15);
% xlabel('estimation error');
% ylabel('prpbability');
% title([testType  ' bswth:' num2str(baseWidth)  ' smwth:' num2str(smoothingWidth) ...
%     '(Hz) Uf:' num2str(boundary.upper) '(Hz) Lf:' num2str(boundary.lower) '(Hz)']);
% outName = ['cumDistFull' num2str(now,30) '.eps'];
% print('-deps',outName);


