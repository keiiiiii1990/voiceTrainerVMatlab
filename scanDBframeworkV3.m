%%  Age DB scan routine
%   by Hideki Kawahara
%   23/May/2012
%   30/June/2012

%dataDirRoot = '/Volumes/Homes-3/Share/Database/Speech_DB/JVPD/';
%dataDirRoot = '//Users/kawahara/Music/JVPD/';
dataDirRoot = '/Volumes/kennkyuuyou_TOSHIBA EXT/0.s130043/database/JVPD/';

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

%%  simple reading loop check


fileCount = 0;
fileCountLimit = 10;
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
                        nameList(fileCount).name = dataName;
                        fullDataPath = [dataDirRoot name '/' dataName];
                        [x,fs] = wavread(fullDataPath); 
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

%%  Test for template generation

outDirectoryName = ['vowels' datestr(now,30)];
mkdir(outDirectoryName);
fileCount = 0;
fileCountLimit = 400;
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
                        nameList(fileCount).name = dataName;
                        fullDataPath = [dataDirRoot name '/' dataName];
                        %disp(fullDataPath);
                        [x,fs] = wavread(fullDataPath);
                        segmentData = vowelSegmentationDB(x,fs);
                        if segmentData.validity
                            disp(dataList(jj).name);
                            vowelTemplate = vowelTemplateGeneratorForDB(segmentData);
                            vowelTemplate.dataName = dataName;
                            save([outDirectoryName '/' dataName(1:end-4)],'vowelTemplate');
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
    if fileCount > fileCountLimit
        break;
    end;
end;
