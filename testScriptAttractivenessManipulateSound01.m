
%%
directoryClass = {'AttractivenessTruncatedVowel2014112816012','UnAttractivenessTruncatedVowel2014112816058'};
directoryCls = {'Attractiveness','UnAttractiveness'};
dirPath = './voice/'; 
for jj = 1:length(directoryClass)
    dncls = directoryClass{jj};
    dcls = directoryCls{jj};
    ClsdirPath = [dirPath dncls];
    voiceDirectoryList = dir(ClsdirPath);
    voiceDirListLen = length(voiceDirectoryList);
    saveExPrmDirName = ['exPrm/' dcls];
    mkdir(saveExPrmDirName);
    for ii = 1:voiceDirListLen
        if strfind(voiceDirectoryList(ii).name,'.wav')
            loadWavFileName = [ClsdirPath '/' voiceDirectoryList(ii).name];
            [dictX,dictFs] = wavread(loadWavFileName);
            extParam = extractPrameterForAttractivenessV01(dictX,dictFs);
            tmpstr = strrep(voiceDirectoryList(ii).name,'.wav','');
            saveExPrmName = [saveExPrmDirName '/' tmpstr 'extPrm.mat'];
            save(saveExPrmName,'extParam');
        end
    end
end

%% —]—Í‚ª‚ ‚ê‚Î‚Â‚­‚é
pairedExPrmDirCls = {'Attractiveness','UnAttractiveness'};
dirPath = './exPrm/';
for ii = 1:length(exprmDirCls)
    dprmCls = exprmDirCls{ii};
    prmClsPth = [dirPath dprmCls '/'];
    prmDirectoryList = dir(prmClsPth);
    
end






