%screeningTest201412

cfgDir = 'config';
cfgDirList = dir(cfgDir);

dirLen = length(cfgDirList);

for ii = 1:dirLen
    if strfind(cfgDirList(ii).name,'.txt')
        %fflg = evaluationVowelSpeech('training',cfgDirList(ii).name);
        evaluationVowelSpeech('training',cfgDirList(ii).name);
    end
end