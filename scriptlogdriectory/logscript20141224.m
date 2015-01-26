%好感度別でもとまった特徴量から変化量をもとめ，音声に適用し，操作後の音声を保存するまでのスクリプトメモ

%20141224version

mkdir 'exPrm02';
mkdir './exPrm02/Attractiveness';

dirpath = '/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141223/exPrm02/';
%attractiveness
[x,fs] = wavread('./voice/AttractivenessTruncatedVowel20141215161714/Attractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV02(x,fs);
save([dirpath 'Attractiveness/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');

%unatt
[x,fs] = wavread('./voice/UnAttractivenessTruncatedVowel20141215161742/UnAttractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV02(x,fs);
save([dirpath 'UnAttractiveness/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');



mkdir 'manipulationDictionary2014122401'

att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141223/exPrm02/Attractiveness/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat');
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141223/exPrm02/UnAttractiveness/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat');
manipulationStructure = manipulationParameterFromDictionaryV02(att,uatt);
save('./manipulationDictionary2014122401/hmdic.mat','manipulationStructure')

%20141215version
%%add 20141223

%create exPrm

mkdir 'exPrm';
mkdir './exPrm/Attractiveness';

dirpath = '/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141223/exPrm02/';
%attractiveness
[x,fs] = wavread('./voice/AttractivenessTruncatedVowel20141215161714/Attractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'Attractiveness/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');

%unatt
[x,fs] = wavread('./voice/UnAttractivenessTruncatedVowel20141215161742/UnAttractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');



mkdir 'manipulationDictionary2014121501'

att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141215/exPrm/Attractiveness/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141215/exPrm/UnAttractiveness/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
manipulationStructure = manipulationParameterFromDictionary(att,uatt);
save('./manipulationDictionary2014121501/hmdic.mat','manipulationStructure')


[x,fs] = wavread('./voice/m1705.wav');
[x,fs] = wavread('./voice/m2801.wav');
%[x,fs] = wavread('./voice/m4201.wav');

clear dictionary
dictionary = load('./manipulationDictionary2014121501/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705DicHM','2014121501m1705eachParameterHM')


dictionary = load('./manipulationDictionary2014121501/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801DicHM','2014121501m2801eachParameterHM')

dictionary = load('./manipulationDictionary2014121501/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201DicHM','2014121501m4201eachParameterHM')

%old version
clear

%create exPrm

mkdir 'exPrm';
mkdir './exPrm/Attractiveness';
dirpath = '/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/';
[x,fs] = wavread('./voice/Attractiveness02TruncatedVowel20141210131848/Attractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');

load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness/AttractivenessAttractiveness_Vowel01_all_TUtrVTTrVextPrm.mat')
att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness/AttractivenessAttractiveness_Vowel01_all_TUtrVTTrVextPrm.mat')
load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness/UnAttractivenessUnAttractiveness_Vowel01_all_TUtrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness/UnAttractivenessUnAttractiveness_Vowel01_all_TUtrVTTrVextPrm.mat')
clear extParam
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
save('./manipulationalDictionary20141202/tudic.mat','manipulationStructure')
clear

[x,fs] = wavread('./voice/m1705.wav');
%[x,fs] = wavread('./voice/m2801.wav');
%[x,fs] = wavread('./voice/m4201.wav');

dictionary = load('./manipulationalDictionary20141202/tsdic.mat')
mSS = AttractiveManipultorByDictionaryVersionLPCSound(x,fs,dictionary);
save('./mSS20141202/m1705dicTS.mat','mSS')
clear dictionary mSS

%audio write
load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/mSS20141202/m1705dicTS.mat')
audiowrite('./modDicSound/m1705dicTSall.wav',mSS.all.outputBuffer,8000);
clear


%20141203
%generateIndivisualCombinationSoundのusage
[x,fs] = wavread('./voice/m1705.wav');
dictionary = load('./manipulationalDictionary20141202/tsdic.mat')
generateIndivisualCombinationSound(x,fs,dictionary,'m1705DicHM','20141203eachParameterHM')


%20141205
[x,fs] = wavread('./voice/m1705.wav');
dictionary = load('./manipulationalDictionary2014120302/hmdic.mat');
generateIndivisualCombinationSound(x,fs,dictionary,'m1705DicHM','2014120501eachParameterHM')

%20141207

%%
dictionary = load('./manipulationalDictionary2014120701/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dichm','201412050702eachParameterHM')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/hsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dichs','201412050702eachParameterHS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/rydic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dicry','201412050702eachParameterRY')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/thdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dicth','201412050702eachParameterTH')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dicts','201412050702eachParameterTS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tudic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m1705Dictu','201412050702eachParameterTU')
%
dictionary = load('./manipulationalDictionary2014120701/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dichm','201412050701m2801eachParameterHM')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/hsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dichs','201412050701m2801eachParameterHS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/rydic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dicry','201412050701m2801eachParameterRY')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/thdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dicth','201412050701m2801eachParameterTH')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dicts','201412050701m2801eachParameterTS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tudic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m2801Dictu','201412050701m2801eachParameterTU')
%
dictionary = load('./manipulationalDictionary2014120701/hmdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dichm','201412050701m4201eachParameterHM')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/hsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dichs','201412050701m4201eachParameterHS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/rydic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dicry','201412050701m4201eachParameterRY')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/thdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dicth','201412050701m4201eachParameterTH')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tsdic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dicts','201412050701m4201eachParameterTS')
clear dictionary
dictionary = load('./manipulationalDictionary2014120701/tudic.mat');
generateIndivisualCombinationSoundV2(x,fs,dictionary,'m4201Dictu','201412050701m4201eachParameterTU')

%20141210
truncatedVowelForAttractiveVowelDB
truncatedVowelForAttractiveVowelDB('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/voice','Attractiveness02')
truncatedVowelForAttractiveVowelDB('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/voice','UnAttractiveness02')
clc
load('AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
extParam
clear
clc
[x,fs] = wavread('./voice/Attractiveness02TruncatedVowel20141210131848/Attractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
exPram
save([dirpath 'Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');
dirpath = '/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/';
save([dirpath 'Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');
clear exPram
clear x fs
[x,fs] = wavread('./voice/Attractiveness02TruncatedVowel20141210131848/Attractiveness_Vowel01_all_HStrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat'],'exPram');
clear x fs dirpath exPram
dirpath = '/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/';
[x,fs] = wavread('./voice/UnAttractiveness02TruncatedVowel20141210131938/UnAttractiveness_Vowel01_all_HMtrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat'],'exPram');
clear x fs exPram
[x,fs] = wavread('./voice/UnAttractiveness02TruncatedVowel20141210131938/UnAttractiveness_Vowel01_all_HStrVTTrV.wav');
exPram = extractPrameterForAttractivenessV01(x,fs);
save([dirpath 'UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat'],'exPram');
clear
load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat')
clear
att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness/UnAttractivenessUnAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
att
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
uatt
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness02/UnAttractivenessUnAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
mkdir 'manipulationalDictionary2014121001'
save('./manipulationalDictionary2014121001/hmdic.mat','manipulationStructure')
att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HMtrVTTrVextPrm.mat')
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
save('./manipulationalDictionary2014121001/hmdic.mat','manipulationStructure')
att = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/Attractiveness02/AttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat')
uatt = load('/Users/amlab/Desktop/StuOnFavVoice/workspace/work20141128/exPrm/UnAttractiveness02/UnAttractivenessAttractiveness_Vowel01_all_HStrVTTrVextPrm.mat')
manipulationStructure = manipulationParameterFromDictionary(att,uatt)
save('./manipulationalDictionary2014121001/hsdic.mat','manipulationStructure')
clc
