%%  Scan template data and play sound
%   by Hideki Kawahara
%   28/Oct./2012

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

%%

templateID = 0;
templateNameList = struct;
sdList = zeros(numberOfTemplate,1);
idList = zeros(numberOfTemplate,1);
for ii = 1:length(templateDirList)
    templateName = templateDirList(ii).name;
    if ~isempty(strfind(templateName,'.mat'))
        templateID = templateID+1;
        templateNameList(templateID).name = templateName;
        tmp = load([dataDirectory '/' templateName]);
        x = tmp.vowelTemplate.signal;
        fs = tmp.vowelTemplate.samplingFrequency;
        disp(templateName)
        for jj = 1:5
            baseIndex = ceil(tmp.vowelTemplate.segmentData.segmentList(jj,1)*fs): ...
                floor(tmp.vowelTemplate.segmentData.segmentList(jj,2)*fs);
%             soundsc(x(max(1,min(length(x),baseIndex))),fs);
        end;
    end;
end;

%%

        tmp = load([dataDirectory '/' templateName]);
        x = tmp.vowelTemplate.signal;
        fs = tmp.vowelTemplate.samplingFrequency;
        disp(templateName)
        for jj = 1:5
            baseIndex = ceil(tmp.vowelTemplate.segmentData.segmentList(jj,1)*fs): ...
                floor(tmp.vowelTemplate.segmentData.segmentList(jj,2)*fs);
%             soundsc(x(max(1,min(length(x),baseIndex))),fs);
        end;

