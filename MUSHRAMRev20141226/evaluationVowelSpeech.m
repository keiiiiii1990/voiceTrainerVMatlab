function evaluationVowelSpeech(varargin)
%mushra for gui Script

%finishFlg

phase=0;
%phasename = 'training';
run_all=false;
random_expe=true;
if length(varargin) > 0,
    phasename = varargin{1};
    if strcmp(phasename,'training'),
        run_all=true;
        if length(varargin) > 1
            cfDirName = varargin{2};
        end
    elseif strcmp(phasename,'evaluation'),
        phase=1;
    end
else
    finishFlg = 1;
end

cfDirName = './config/mushram_configALL.txt';

%if strfind(cfgDirList(ii).name,'.txt')
%fid=fopen('mushram_config.txt','r');
%fid=fopen(cfgDirList(ii).name,'r');
%fid=fopen('mushram_config.txt','r');
%fid=fopen('mushram_configMALE.txt','r');
%fid=fopen('mushram_configFEMALE.txt','r');
fid=fopen(cfDirName,'r');
config=fscanf(fid,'%c');
fclose(fid);

if length(config)>1,
    config=strrep(config,char(13),char(10));                                %char(13):CR(windows) to char(10):LF(unix), 
    c=find(config~=10);
    config=[config(1:c(end)) char(10) char(10)];
    while ~isempty(strfind(config,[char(10) char(10) char(10)])),
        config=strrep(config,[char(10) char(10) char(10)],[char(10) char(10)]);
    end
else,
    errordlg('The configuration file mushram_config.txt is empty.','Error');
    return;
end



dblines=strfind(config,[char(10) char(10)]);                                %改行2回の位置
nbexpe=length(dblines);
expconfig=config(1:dblines(1));
lines=strfind(expconfig,char(10));                                          %実験一回分のデータ数（改行位置）


%files=cell(nbexpe,nbfile);
dblines=[-1 dblines];

%{
if exist
if handles.expe > 2 
    expconfig = config(dblines(handles.expe)+2:dblines(handles.expe+1));
    lines = strfind(expconfig,char(10));
    nbfile=length(lines);
end
end
%}

%%{
%ファイル数の多い実験のファイル数を求める．
larLinVal = 0;
fileCntEchExp = zeros(1,nbexpe);

teConf = expconfig;
for ii = 1:length(dblines)-1
    teConf = config(dblines(ii)+2:dblines(ii+1));
    fileCntEchExp(ii) = length(strfind(teConf,char(10)));
    if larLinVal < length(strfind(teConf,char(10)))
        larLinVal = length(strfind(teConf,char(10)));
    end
end
nbfile = larLinVal;
files=cell(nbexpe,larLinVal);

%}

%よくわからない処理
for e=1:length(dblines)-1,                                                  
    expconfig=config(dblines(e)+2:dblines(e+1));                            %1回分の実験データ記録, 最初は2文字分（改行分?）飛ばされている．
    lines=strfind(expconfig,char(10));
    %if length(lines) == nbfile,
    lines=[0 lines];
    %for f=1:length(lines)-1,
    for f = 1:fileCntEchExp(e)
        file=expconfig(lines(f)+1:lines(f+1)-1);
        if exist(file),
            files{e,f}=file;
        else,
            errordlg(['The specified sound file ' file ' does not exist. Check the configuration file mushram_config.txt.'],'Error');
            return;
        end
    end
    %{    
    else,
        errordlg('The number of test files must be the same for all experiments. Check the configuration file mushram_config.txt.','Error');
        return;
    end
    %}
end

%たぶんランダム表示
expe_order=randperm(nbexpe);
file_order=randperm(nbfile);

%%%randomizing the order of the experiments or checking the alternative experiment order is acceptable
if exist('expe_order'),
    err=false;
    for e=1:nbexpe,
        if ~any(expe_order==e),
            err=true;
        end
    end
    if err,
        errordlg('Bad input parameters.','Error');
        return;
    end
elseif random_expe,
    expe_order=randperm(nbexpe);
else,
    expe_order=1:nbexpe;    
end

if phase,
    resSaveFile = 'mushram_results';
    [filename,pathname]=uiputfile([resSaveFile datestr(now,30) '.txt'],'Results file name');
    resultfile=[pathname filename];

    if ~resultfile,
        return;
    end

    fig=evaluation_gui(nbfile,1,nbexpe);
    if ~fig,
        errordlg('There are too many test files to display. Try increasing the resolution of the screen.','Error');
        return;
    end
    handles=guihandles(fig);

    %randomizing the order of the tested files
    file_order=randperm(nbfile);

    %storing data within the GUI
    handles.expe=1;
    handles.ratings=zeros(nbexpe,nbfile);
    handles.resultfile=resultfile;
else,
    %fig=training_gui(nbfile,nbexpe,run_all);
    fig=training_gui(fileCntEchExp,nbexpe,run_all);

    if ~fig,
        errordlg('There are too many experiments or too many test files to display. Try increasing the resolution of the screen.','Error');
        return;
    end

    handles=guihandles(fig);

    %randomizing the order of the tested files
    
    %{
    file_order=zeros(nbexpe,nbfile-1);
    for e=1:nbexpe,
        file_order(e,:)=randperm(nbfile-1)+1;
    end
    %}
    
    file_order=zeros(nbexpe,nbfile);
    for e=1:nbexpe,
        file_order(e,:)=randperm(nbfile);
    end
    
    %file_order=randperm(nbfile);
    %storing data within the GUI
    handles.run_all=run_all;
end

handles.expe_order=expe_order;
handles.nbexpe=nbexpe;
handles.file_order=file_order;
handles.nbfile=nbfile;
handles.files=files;
handles.fileCntEchExp = fileCntEchExp;
handles.time=clock;
guidata(fig,handles);



