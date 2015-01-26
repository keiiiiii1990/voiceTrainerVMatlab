function truncatedVowelForAttractiveVowelDB(OriginalVoiceDirName,Opt)
OVDir = OriginalVoiceDirName;
Cls = Opt;% attractivenss or Unattractiveness

%% make vowel template memoScript
%/Users/amlab/MATLABFORSTRAIGHT/mfiles/で実行が必要．

%これらのファイルを/Users/amlab/MATLABFORSTRAIGHT/mfiles/へ

saveVSDirName = ['VowelStucture' num2str(fix(clock))];
mkdir(saveVSDirName);
saveVTDirName = ['VowelTemplate' num2str(fix(clock))];
mkdir(saveVTDirName)

OVDirAddClass = [OVDir '/' Opt]
dirList = dir(OVDirAddClass)
dirListLen = length(dirList);

for ii = 1:dirListLen
    if strfind(dirList(ii).name,'.wav')
        
        %generate Vowel Structure
        saveVSName = strrep(dirList(ii).name,'.wav','');
        LoadWavFileName = [OVDirAddClass '/' dirList(ii).name];
        [x,fs] = wavread(LoadWavFileName);
        vs = vowelSegmentationSimpleV3(x,fs);
        saveVSFileName = [saveVSDirName '/' saveVSName 'VS.mat'];
        save(saveVSFileName,'vs');
        
        %make vowelTemplate
        segmentData = vs;
        vowelTemplate = vowelTemplateGeneratorForDBV3(segmentData);
        vowelTemplate.dataName = dirList(ii).name;
        vowelTemplate.signal = x;
        vowelTemplate.samplingFrequency = fs;
        vowelTemplate.gender = 'm';
        
        saveVTName = strrep(dirList(ii).name,'.wav','');
        saveVTFileName = [saveVTDirName '/' saveVTName 'VT.mat'];
        save(saveVTFileName,'vowelTemplate')
    end
end


%% trucate vowel segmentation
%ref:memoOnTemplate-4.pdf
%ref:/Users/amlab/MATLABFORSTRAIGHT/mfiles/vowelSegmentationSimpleV3.m

%load vowel template MAT file

VTDirList = dir(saveVTDirName);
VTDirListLen = length(VTDirList);

saveTrWavDir = [Cls 'TruncatedVowel' num2str(fix(clock))];
mkdir(saveTrWavDir);

for jj = 1:VTDirListLen
    if strfind(VTDirList(jj).name,'.mat')
        LoadMatFileName = [saveVTDirName '/' VTDirList(jj).name];
        value = load(LoadMatFileName);
        repMatName = strrep(VTDirList(jj).name,'.mat',''); 
        saveWavFileName = [saveTrWavDir '/' repMatName 'TrV.wav'];
        fs = value.vowelTemplate.segmentData.samplingFrequency;
        y = value.vowelTemplate.segmentData.waveform;
        vowelSpecInfo = struct;
        for ii = 1:5
            % choice segment range
            yStart = value.vowelTemplate.segmentData.segmentList(ii, 1);
            yEnd = value.vowelTemplate.segmentData.segmentList(ii, 2);
            yUse = y(round(yStart*fs):round(yEnd*fs));
            yUse2 = yUse.*hamming(length(yUse));
            yUse3 = zeros(1, length(yUse2));
            yUse3(1) = 0;
            for l = 2:length(yUse3)
                yUse3(l) = yUse(l) - 0.98*yUse(l-1); % preEmphasis
            end;
            switch ii
                case 1
                    vowelSpecInfo.original.waveform.a = yUse2;
                case 2
                    vowelSpecInfo.original.waveform.i = yUse2;
                case 3
                    vowelSpecInfo.original.waveform.u = yUse2;
                case 4
                    vowelSpecInfo.original.waveform.e = yUse;
                case 5
                    vowelSpecInfo.original.waveform.o = yUse2;
            end;
        end;
        vowelSpecInfo.original.samplingFrequency = fs;

        alength = length(vowelSpecInfo.original.waveform.a);
        ilength = length(vowelSpecInfo.original.waveform.i);
        ulength = length(vowelSpecInfo.original.waveform.u);
        elength = length(vowelSpecInfo.original.waveform.e);
        olength = length(vowelSpecInfo.original.waveform.o);
        trvowelWaveform = zeros(alength+ilength+ulength+elength+olength,1);
        nextend = alength;
        trvowelWaveform(1:alength) = vowelSpecInfo.original.waveform.a;

        prenextend = nextend;
        nextend = nextend + ilength;
        trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.i;

        prenextend = nextend;
        nextend = nextend + ulength;
        trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.u;

        prenextend = nextend;
        nextend = nextend + elength;
        trvowelWaveform(prenextend+1:nextend) = vowelSpecInfo.original.waveform.e;

        prenextend = nextend;
        nextend = nextend + olength;
        trvowelWaveform(prenextend+1:end) = vowelSpecInfo.original.waveform.o;
        audiowrite(saveWavFileName,trvowelWaveform,vowelSpecInfo.original.samplingFrequency);
    end
end
        

end