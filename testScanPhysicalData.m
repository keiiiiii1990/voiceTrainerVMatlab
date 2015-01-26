%%  Physical data check
%   by Hideki Kawahara
%   04/Nov./2012


dataDirectory = 'vowels20121028T185523OK';

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

load smoothingTest1

dbOut = scanPhysicalData('rawPhysicalinfoEdt.txt');

%%  check name agreement

matchCount = 0;
indexList = zeros(length(dbOut.height),1);
for ii = 1:length(dbOut.height)
    templateName = char(dbOut.name{ii});
    for jj = 1:length(dbOut.height)
        dbName = templateNameList(jj).name;
        if strcmp(templateName,dbName(1:5))
            indexList(ii) = jj;
            matchCount = matchCount+1;
            break;
        end;
    end;
end;
if matchCount == length(dbOut.height)
    disp('names are all arigned');
end;

alignedVTL = dbOut.height*0;
for ii = 1:length(dbOut.height)
    alignedVTL(ii) = distanceStructureBase.vtl(indexList(ii));
end;

%%

figure;
femaleHightSrt = sort(dbOut.height(dbOut.height~=0 & dbOut.gender==0));
maleHightSrt = sort(dbOut.height(dbOut.height~=0 & dbOut.gender==1));
plot(femaleHightSrt,(1:length(femaleHightSrt))/length(femaleHightSrt),'r','linewidth',2);
hold on;
plot(maleHightSrt,(1:length(maleHightSrt))/length(maleHightSrt),'b','linewidth',2);
grid on;
set(gca,'fontsize',15);
xlabel('hight');
ylabel('probability');
legend('female','male');
outName = ['hightdist' datestr(now,30) '.eps'];
print('-depsc',outName);

%%

figure;
femaleWeightSrt = sort(dbOut.weight(dbOut.weight~=0 & dbOut.gender==0));
maleWeightSrt = sort(dbOut.weight(dbOut.weight~=0 & dbOut.gender==1));
plot(femaleWeightSrt,(1:length(femaleWeightSrt))/length(femaleWeightSrt),'r','linewidth',2);
hold on;
plot(maleWeightSrt,(1:length(maleWeightSrt))/length(maleWeightSrt),'b','linewidth',2);
grid on;
set(gca,'fontsize',15);
xlabel('weight');
ylabel('probability');
legend('female','male');
outName = ['weightdist' datestr(now,30) '.eps'];
print('-depsc',outName);

%%
figure;
femaleVTLSrt = sort(alignedVTL(dbOut.height~=0 & dbOut.gender==0));
maleVTLSrt = sort(alignedVTL(dbOut.height~=0 & dbOut.gender==1));
plot(femaleVTLSrt,(1:length(femaleVTLSrt))/length(femaleVTLSrt),'r','linewidth',2);
hold on;
plot(maleVTLSrt,(1:length(maleVTLSrt))/length(maleVTLSrt),'b','linewidth',2);
grid on;
set(gca,'fontsize',15);
xlabel('relative VTL');
ylabel('probability');
legend('female','male');
outName = ['VTLdistN' datestr(now,30) '.eps'];
print('-depsc',outName);

%%

femaleHightUSrt = dbOut.height(dbOut.height~=0 & dbOut.gender==0);
maleHightUSrt = dbOut.height(dbOut.height~=0 & dbOut.gender==1);
femaleVTLUSrt = alignedVTL(dbOut.height~=0 & dbOut.gender==0);
maleVTLUSrt = alignedVTL(dbOut.height~=0 & dbOut.gender==1);
figure;
plot(femaleVTLUSrt,femaleHightUSrt,'ro');
hold on;
plot(maleVTLUSrt,maleHightUSrt,'bo');
grid on;
lsline;
set(gca,'fontsize',15);
ylabel('height');
xlabel('relative VTL');
legend('female','male');
outName = ['VTLvsHeightSctr' datestr(now,30) '.eps'];
print('-depsc',outName);

% %%
figure;
femaleWeightUSrt = (dbOut.weight(dbOut.weight~=0 & dbOut.gender==0));
maleWeightUSrt = (dbOut.weight(dbOut.weight~=0 & dbOut.gender==1));

plot(femaleVTLUSrt,femaleWeightUSrt,'ro');
hold on;
plot(maleVTLUSrt,maleWeightUSrt,'bo');
grid on;
lsline;
set(gca,'fontsize',15);
ylabel('weight');
xlabel('relative VTL');
legend('female','male');
outName = ['VTLvsWeightSctr' datestr(now,30) '.eps'];
print('-depsc',outName);


