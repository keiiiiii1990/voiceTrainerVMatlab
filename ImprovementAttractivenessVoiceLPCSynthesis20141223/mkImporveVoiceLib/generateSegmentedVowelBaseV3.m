function vowelBaseStatus = generateSegmentedVowelBaseV3(dataDirRoot,fileCountLimit)
%   vowelBaseStatus = generateSegmentedVowelBaseV2(dataDirRoot,fileCountLimit)
%       Input
%           dataDirRoot     : JVPD database location
%           fileCountLimit  : for debug

%   Designed and coded by Hideki Kawahara
%   28/Oct./2012 extended segment range with power and f0 : V2
%   28/Oct./2012 extended segment range with power, f0 and spDstnc : V3

dataDirRoot = '//media.sys.wakayama-u.ac.jp/s130043/s130043/M1/JVPD/';

dirContents = dir(dataDirRoot);

fileCount = 0;
for ii = 1:length(dirContents)
    name = dirContents(ii).name;
    switch name(1)
        case {'.' '..'}
        case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
            %disp(dirContents(ii).name);
            dataList = dir([dataDirRoot dirContents(ii).name]);
            for jj = 1:length(dataList)
                dataName = dataList(jj).name;
                switch dataName(1)
                    case {'f' 'm'}
                        %disp(dataList(jj).name);
                        fileCount = fileCount+1;
                end;
            end;
    end;
end;
numberOfFiles = fileCount;

%%  Test for template generation

outDirectoryName = ['vowels' datestr(now,30)];
mkdir(outDirectoryName);
%fileCount = 0;
fileCountLimit = numberOfFiles; % just for test
vowelBaseStatus = struct;
for ii = 1:length(dirContents)
    name = dirContents(ii).name;
    switch name(1)
        case {'.' '..'}
        case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
            %disp(dirContents(ii).name);
            dataList = dir([dataDirRoot dirContents(ii).name]);
            for jj = 1:length(dataList)
                dataName = dataList(jj).name;
                switch dataName(1)
                    case {'f' 'm'}
                        %disp(dataList(jj).name);
                        fileCount = fileCount+1;
                        %nameList(fileCount).name = dataName;
                        fullDataPath = [dataDirRoot name '/' dataName];
                        [x,fs] = wavread(fullDataPath); 
                        segmentData = vowelSegmentationSimpleV3(x,fs);
                        if segmentData.validity
                            disp(dataList(jj).name);
                            vowelTemplate = vowelTemplateGeneratorForDBV3(segmentData);
                            vowelTemplate.dataName = dataName;
                            vowelTemplate.signal = x;
                            vowelTemplate.samplingFrequency = fs;
                            vowelTemplate.gender = dataName(1);
                            vowelTemplate.age = str2double(char(dataName(2:3)));
                            vowelTemplate.ID = str2double(char(dataName(4:5)));
                            save([outDirectoryName '/' dataName(1:end-4)],'vowelTemplate');
                        end;
                        vowelBaseStatus.fileName{fileCount} = dataName;
                        vowelBaseStatus.validity{fileCount} = segmentData.validity;
                        if fileCount > fileCountLimit
                            break;
                        end;
                end;
                if fileCount > fileCountLimit
                    break;
                end;
            end;
            if fileCount > fileCountLimit
                break;
            end;
    end;
    if fileCount > fileCountLimit
        break;
    end;
end;
%%
%vowelBaseStatus.nameList = nameList;
return;
