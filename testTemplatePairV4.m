%%  Test sctript for vowel template pairing
%   Originally designed by Hideki Kawahara
%   24/May/2012
%   25/May/2012 test elimination outlier
%   30/June/2012 test all data
%   14/July/2012 bind annotation and analyisis results

dataDirectory = 'vowels20121114T221850';
templateDirList = dir(dataDirectory);

%%  generate template name list

templateID = 0;
templateNameList = struct;
for ii = 1:length(templateDirList)
    templateName = templateDirList(ii).name;
    if ~isempty(strfind(templateName,'.mat'))
        templateID = templateID+1;
        templateNameList(templateID).name = templateName;
    end;
end;
numberOfTemplate = templateID;

%%  generate 3 dimentional raw template

tmp = load([dataDirectory '/' templateNameList(1).name]);
frequencyAxis = tmp.vowelTemplate.frequencyAxis;
rawSize = length(frequencyAxis);
rawTemplate = zeros(rawSize,5,numberOfTemplate);
rawF0List = zeros(5,numberOfTemplate);

for ii = 1:numberOfTemplate
    tmp = load([dataDirectory '/' templateNameList(ii).name]);
    rawTemplate(:,:,ii) = tmp.vowelTemplate.rawTemplate;
    rawF0List(:,ii) = tmp.vowelTemplate.segmentMedianF0;
end;

%%  check for exsisting calculation

if ~exist('snap20120716.mat')
    
    %%  template distances
    
    baseWidth = 2500;
    smoothingWidth = 300;
    boundary.upper = 5000;
    boundary.lower = 300;
    
    distanceMatrix = zeros(numberOfTemplate,numberOfTemplate);
    VTLRatioMatrix = zeros(numberOfTemplate,numberOfTemplate);
    temlateStrA.frequencyAxis = frequencyAxis(:);
    temlateStrB.frequencyAxis = frequencyAxis(:);
    for ii = 1:numberOfTemplate
        disp([num2str(ii) ' ' datestr(now)]);
        temlateStrA.vowel = rawTemplate(:,:,ii);
        temlateStrA.F0 = rawF0List(:,ii);
        for jj = 1:numberOfTemplate
            temlateStrB.vowel = rawTemplate(:,:,jj);
            temlateStrB.F0 = rawF0List(:,jj);
            tr = compareTemplatesDB(temlateStrA,temlateStrB,baseWidth,smoothingWidth,boundary);
            distanceMatrix(ii,jj) = tr.minvalue;
            VTLRatioMatrix(ii,jj) = tr.bestRatio;
        end;
    end;
    save('snap20120716')
else
    load snap20120716.mat
end;

%%
nameCell = cell(numberOfTemplate);
for ii = 1:numberOfTemplate
    nameCell{ii} = templateNameList(ii).name(1:end-4);
end;
%%  print results

figure(1);imagesc(VTLRatioMatrix);colorbar
%set(gca,'YTickLabel',nameCell,'fontsize',13)
set(gca,'fontsize',15)
xlabel('tentative talker ID');
ylabel('tentative talker ID');
title(['VTL ratio: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(1),'VTLratio.fig');

figure(2);imagesc(distanceMatrix);colorbar
%set(gca,'YTickLabel',nameCell,'fontsize',13)
set(gca,'fontsize',15)
xlabel('tentative talker ID');
ylabel('tentative talker ID');
title(['distance matrix: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(2),'DistanceRatio.fig');

figure(3);imagesc(abs(log(VTLRatioMatrix)));colorbar
%set(gca,'YTickLabel',nameCell,'fontsize',13)
set(gca,'fontsize',15)
xlabel('tentative talker ID');
ylabel('tentative talker ID');
title(['abs of log(VTL ratio): ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(3),'AbsOfLog_VTLratio.fig');

figure(4);plot(VTLRatioMatrix(:),distanceMatrix(:),'.','markersize',12);grid on;
set(gca,'fontsize',15)
xlabel('VTL ratio');
ylabel('minimum distance (dB)');
title(['scatter plot: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(4),'ScatterPlot.fig');

%%

H = zeros(numberOfTemplate*(numberOfTemplate-1)+1,numberOfTemplate);
V = zeros(numberOfTemplate*(numberOfTemplate-1)+1,1);
rowID = 0;
for ii = 1:numberOfTemplate
    for jj = 1:numberOfTemplate
        if ii ~= jj
            rowID = rowID+1;
            H(rowID,ii) = 1;
            H(rowID,jj) = -1;
            V(rowID) = -log(VTLRatioMatrix(ii,jj));
        end;
    end;
end;
rowID = rowID+1;
H(rowID,:) = 1;
V(rowID) = 0;
VTL = inv(H'*H)*(H'*V);
vtl = exp(VTL);

estimatedVTLR = zeros(numberOfTemplate,numberOfTemplate);
for ii = 1:numberOfTemplate
    for jj = 1:numberOfTemplate
        estimatedVTLR(ii,jj) = vtl(jj)/vtl(ii);
    end;
end;

%%

figure(5);plot(estimatedVTLR(:),VTLRatioMatrix(:),'.','markersize',12);grid on;
set(gca,'fontsize',15)
xlabel('VTL ratio by regression analysis');
ylabel('spectral VTL ratio');
title(['regression: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(5),'Regression.fig');

%%
[valuev,sortedIndexv] = sort(vtl);
estimationErrorList = zeros(numberOfTemplate,1);
for ii = 1:numberOfTemplate
    tmpError = 0;
    for jj = 1:numberOfTemplate
        tmpError = tmpError+(estimatedVTLR(ii,jj)-VTLRatioMatrix(ii,jj)).^2;
    end;
    tmpError = sqrt(tmpError/numberOfTemplate);
    estimationErrorList(ii) = tmpError;
    disp(['name:' nameCell{sortedIndexv(ii)} '  VTL:' num2str(vtl(sortedIndexv(ii))) ...
        '  error:' num2str(tmpError) ]);
end;

figure(6);plot(sort(estimationErrorList));grid on;
title(['estimation error: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
ylabel('rms error in VTL ratio');
xlabel('cumulative count');
saveas(figure(6),'EstimationError.fig');

%%
rmsDistanceList = sqrt(mean(distanceMatrix.^2));
%%
figure(7);stem(rmsDistanceList);grid on;
set(gca,'fontsize',14);
xlabel('sorted talkers re. error');
ylabel('rms error in VTL re. regression VTL');
title('regression VTL to spectral VTL');
saveas(figure(7),'RegressionVTLtoSpectralVTL.fig');

figure(8);plot(sort(rmsDistanceList));grid on;
set(gca,'fontsize',14);
xlabel('sorted talkers re. distance');
ylabel('rms spectral distance');
title('minimum spectral distance to all talkers');
saveas(figure(8),'MinimumSpectralDistance.fig');

%%  repeat analysis on reliable data only

%spectral distances smaller than 8dB can be reliable data
numberOfReliableOnes = sum(rmsDistanceList<8);
indexList = 1:numberOfTemplate;
reliableIndex = indexList(rmsDistanceList<8);

H = zeros(numberOfReliableOnes*(numberOfReliableOnes-1)+1,numberOfReliableOnes);
V = zeros(numberOfReliableOnes*(numberOfReliableOnes-1)+1,1);
rowID = 0;
for ii = 1:numberOfReliableOnes
    for jj = 1:numberOfReliableOnes
        if ii ~= jj
            rowID = rowID+1;
            H(rowID,ii) = 1;
            H(rowID,jj) = -1;
            V(rowID) = -log(VTLRatioMatrix(reliableIndex(ii),reliableIndex(jj)));
        end;
    end;
end;
rowID = rowID+1;
H(rowID,:) = 1;
V(rowID) = 0;
VTL = inv(H'*H)*(H'*V);
vtl = exp(VTL);

estimatedVTLR = zeros(numberOfReliableOnes,numberOfReliableOnes);
for ii = 1:numberOfReliableOnes
    for jj = 1:numberOfReliableOnes
        estimatedVTLR(ii,jj) = vtl(jj)/vtl(ii);
    end;
end;

reliableVTLRatioMatrix = VTLRatioMatrix(rmsDistanceList<8,:);
reliableVTLRatioMatrix = reliableVTLRatioMatrix(:,rmsDistanceList<8);

%%
figure(9);plot(estimatedVTLR(:),reliableVTLRatioMatrix(:),'.','markersize',12);grid on;
set(gca,'fontsize',15)
xlabel('reliable VTL ratio by regression analysis');
ylabel('spectral VTL ratio');
title(['regression: ' num2str(boundary.lower) ' to ' num2str(boundary.upper) ' (Hz)']);
saveas(figure(9),'Regression_reliable.fig');

%%  Read physical data and check for basic statistics

ageList = zeros(numberOfReliableOnes,1);
genderList = zeros(numberOfReliableOnes,1);
for ii = 1:numberOfReliableOnes
    currentName = nameCell{reliableIndex(ii)};
    ageList(ii) = str2num(currentName(2:3));
    switch currentName(1)
        case 'f'
            genderList(ii) = 0;
        otherwise
            genderList(ii) = 1;
    end;
    %genderList(ii) = currentName(1);
end;

%%  display age, gender and vtl regression results

figure(10);plot(ageList(genderList==0),vtl(genderList==0),'.r','markersize',12);grid on;
hold on
plot(ageList(genderList==1),vtl(genderList==1),'.b','markersize',12);grid on;
set(gca,'fontsize',14);
xlabel('age');
ylabel('estimated relative VTL');
title('regression results: raw data');
saveas(figure(10),'RegressionResults_rawdata.fig');

%%  read physical data

dbOut = scanPhysicalData('rawPhysicalinfoEdt.txt');

%%  diaplay physical data

figure(11);plot(dbOut.age(dbOut.gender==0 & dbOut.height>0),...
    dbOut.height(dbOut.gender==0 & dbOut.height>0),'.r','markersize',12);grid on;
hold on;
plot(dbOut.age(dbOut.gender==1 & dbOut.height>0), ...
    dbOut.height(dbOut.gender==1 & dbOut.height>0),'.b','markersize',12);grid on;
set(gca,'fontsize',14);
xlabel('age');
ylabel('height (cm)');
title('physical data');
saveas(figure(11),'PhysicalData1.fig');

figure(12);plot(dbOut.age(dbOut.gender==0 & dbOut.height>0),...
    dbOut.weight(dbOut.gender==0 & dbOut.height>0),'.r','markersize',12);grid on;
hold on;
plot(dbOut.age(dbOut.gender==1 & dbOut.height>0), ...
    dbOut.weight(dbOut.gender==1 & dbOut.height>0),'.b','markersize',12);grid on;
set(gca,'fontsize',14);
xlabel('age');
ylabel('weight (kg)');
title('physical data');
saveas(figure(12),'PhysicalData2.fig');

figure(13);semilogy(dbOut.height(dbOut.gender==0 & dbOut.height>0),...
    dbOut.weight(dbOut.gender==0 & dbOut.height>0),'.r','markersize',12);grid on;
hold on;
semilogy(dbOut.height(dbOut.gender==1 & dbOut.height>0), ...
    dbOut.weight(dbOut.gender==1 & dbOut.height>0),'.b','markersize',12);grid on;
semilogy(100:200,10*((100:200)/100).^2,'--');
semilogy(100:200,10*((100:200)/100).^3,'--');
axis([100 200 10 140]);
set(gca,'fontsize',14);
xlabel('height (cm)');
ylabel('weight (kg)');
title('physical data');
saveas(figure(13),'PhysicalData3.fig');
%%
tic
matchedIndex = zeros(dbOut.numberOfItems,1);
matchedTemplateIndex = zeros(size(templateNameList,2),1);
templateItems = size(templateNameList,2);
for ii = 1:dbOut.numberOfItems
    keyName = char(dbOut.name{ii});
    for jj = 1:templateItems
        testName = char(templateNameList(jj).name);
        if strcmp(keyName,testName(1:5))
            matchedIndex(ii) = jj;
            matchedTemplateIndex(jj) = ii;
        end;
    end;
end;
toc

%%

figure(14);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    switch dbOut.gender(jj)
        case 1
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'b.','markersize',12);end;
        case 0
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for all ages');
saveas(figure(14),'VTL-HeightRegressionForAllAges.fig');

%%

figure(15);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    if dbOut.age(jj) < 20
    switch dbOut.gender(jj)
        case 1
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'b.','markersize',12);end;
        case 0
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
    end;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for young');
saveas(figure(15),'VTL-HeightRegressionForYoung.fig');

%%

figure(16);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    if dbOut.age(jj) >= 20
    switch dbOut.gender(jj)
        case 1
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'b.','markersize',12);end;
        case 0
            if dbOut.height(jj)>0 ;plot(vtl(ii),dbOut.height(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
    end;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for adults');
saveas(figure(16),'VTL-HeightRegressionForAdults.fig');

%%

figure(17);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    switch dbOut.gender(jj)
        case 1
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'b.','markersize',12);end;
        case 0
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for all ages');
saveas(figure(17),'VTL-WeightRegressionForAllAges.fig');

%%

figure(18);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    if dbOut.age(jj) < 20
    switch dbOut.gender(jj)
        case 1
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'b.','markersize',12);end;
        case 0
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
    end;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for young');
saveas(figure(18),'VTL-WeightRegressionForYoung.fig');

%%

figure(19);
for ii = 1:length(reliableIndex)
    jj = matchedTemplateIndex(reliableIndex(ii));
    if dbOut.age(jj) >= 20
    switch dbOut.gender(jj)
        case 1
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'b.','markersize',12);end;
        case 0
            if dbOut.weight(jj)>0 ;plot(vtl(ii),dbOut.weight(jj),'r.','markersize',12);end;
    end;
    grid on;hold on;
    end;
end;
set(gca,'fontsize',14);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for adults');
saveas(figure(19),'VTL-WeightRegressionForAdults.fig');
