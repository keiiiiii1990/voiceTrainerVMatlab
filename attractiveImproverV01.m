function varargout = attractiveImproverV01(varargin)
% ATTRACTIVEIMPROVERV01 MATLAB code for attractiveImproverV01.fig
%      ATTRACTIVEIMPROVERV01, by itself, creates a new ATTRACTIVEIMPROVERV01 or raises the existing
%      singleton*.
%
%      H = ATTRACTIVEIMPROVERV01 returns the handle to a new ATTRACTIVEIMPROVERV01 or the handle to
%      the existing singleton*.
%
%      ATTRACTIVEIMPROVERV01('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ATTRACTIVEIMPROVERV01.M with the given input arguments.
%
%      ATTRACTIVEIMPROVERV01('Property','Value',...) creates a new ATTRACTIVEIMPROVERV01 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before attractiveImproverV01_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to attractiveImproverV01_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help attractiveImproverV01

% Last Modified by GUIDE v2.5 03-Jan-2015 17:47:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @attractiveImproverV01_OpeningFcn, ...
                   'gui_OutputFcn',  @attractiveImproverV01_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before attractiveImproverV01 is made visible.
function attractiveImproverV01_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to attractiveImproverV01 (see VARARGIN)

% Choose default command line output for attractiveImproverV01
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%% mycode
initializeDisplay(handles);

% UIWAIT makes attractiveImproverV01 wait for user response (see UIRESUME)
% uiwait(handles.AttractivenessImprovementMainGUI);

%%%%%%%%%%% my code
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

set(myGUIdata.Start_Recording, 'enable', 'on');

set(myGUIdata.Stop_Recording, 'enable', 'off');
set(myGUIdata.Save_RecordingSound, 'enable', 'off');
set(myGUIdata.Play_modified_Sound, 'enable', 'off');
set(myGUIdata.Play_Original_Sound, 'enable', 'off');
set(myGUIdata.Save_Modified_Sound, 'enable', 'off');
set(myGUIdata.make_Improve_sound, 'enable', 'off');

myGUIdata.smallviewerWidth = 30;
myGUIdata.samplingFrequency = 44100;
myGUIdata.recordObj1 = audiorecorder(myGUIdata.samplingFrequency,24,1);
%record(myGUIdata.recordObj1)
switch get(myGUIdata.recordObj1,'Running')
    case 'on'
        stop(myGUIdata.recordObj1);
end
guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
%Start_Recording_Callback(hObject, eventdata, handles)
end


% --- Outputs from this function are returned to the command line.
function varargout = attractiveImproverV01_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function initializeDisplay(handles)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);
myGUIdata.maxAudioRecorderCount = 200;
myGUIdata.audioRecorderCount = myGUIdata.maxAudioRecorderCount;
myGUIdata.maxLevelIndicator = -100*ones(myGUIdata.maxAudioRecorderCount,1);
myGUIdata.yMax = 1;

fs = 44100;
dataLength = round(30/1000*fs);
fftl = 2.0.^ceil(log2(dataLength));
fAxis = (0:fftl-1)/fftl*fs;
w = blackman(dataLength);
pw = 20*log10(abs(fft(randn(dataLength,1).*w,fftl)/sqrt(sum(w.^2))));


%for LPC processing
myGUIdata.ForLPC.samplingFrequency = 8000; % in Hz
myGUIdata.ForLPC.windowLength = 0.08; % in second
myGUIdata.ForLPC.windowLengthInSamples = round(myGUIdata.ForLPC.windowLength*myGUIdata.ForLPC.samplingFrequency/2)*2+1;
myGUIdata.ForLpc.fftl = 2.0^ceil(log2(myGUIdata.ForLPC.windowLengthInSamples));

%tmp is window
tmp = blackman(myGUIdata.ForLPC.windowLengthInSamples);
myGUIdata.ForLPC.window = tmp/sqrt(sum(tmp.^2));

axes(handles.orgVTAFView);
crossSection =  [0.2803; ...
    0.6663; ...
    0.5118; ...
    0.3167; ...
    0.1759; ...
    0.1534; ...
    0.1565; ...
    0.1519; ...
    0.0878; ...
    0.0737];
[X,Y,Z] = cylinder(crossSection,40);
myGUIdata.tract3DOrg = surf(Z,Y,X);
view(-26,12);
axis([0 1 -1 1 -1 1]);
axis off;
axis('vis3d');
rotate3d on;

axes(handles.modVTAFView);
crossSection =  [0.2803; ...
    0.6663; ...
    0.5118; ...
    0.3167; ...
    0.1759; ...
    0.1534; ...
    0.1565; ...
    0.1519; ...
    0.0878; ...
    0.0737];
[X,Y,Z] = cylinder(crossSection,40);
myGUIdata.tract3DMod = surf(Z,Y,X);
view(-26,12);
axis([0 1 -1 1 -1 1]);
axis off;
axis('vis3d');
rotate3d on;

guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end


% --- Executes on button press in Start_Recording.
function Start_Recording_Callback(hObject, eventdata, handles)
% hObject    handle to Start_Recording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

%off button
set(myGUIdata.Save_RecordingSound, 'enable', 'off');
set(myGUIdata.Start_Recording, 'enable', 'off');
set(myGUIdata.make_Improve_sound, 'enable', 'off');

set(myGUIdata.Save_RecordingSound, 'enable', 'off');
set(myGUIdata.Save_Modified_Sound, 'enable', 'off');
set(myGUIdata.Play_modified_Sound, 'enable', 'off');
set(myGUIdata.Play_Original_Sound, 'enable', 'off');

%on button
set(myGUIdata.Stop_Recording, 'enable', 'on');

%set(myGUIdata.vtlDisplay,'visible','off');
%set(myGUIdata.vtlEstimateButton,'enable','off');
%{
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms);
end
%}

myGUIdata.audioRecorderCount = myGUIdata.maxAudioRecorderCount;
myGUIdata.maxLevelIndicator = -100*ones(myGUIdata.maxAudioRecorderCount,1);
myGUIdata.yMax = 1;
record(myGUIdata.recordObj1);
%{
switch get(myGUIdata.timer50ms,'running')
    case 'off'
        start(myGUIdata.timer50ms);
    case 'on'
    otherwise
        disp('timer is bloken!');
end
%}
guidata(handles.AttractivenessImprovementMainGUI);

end

% --- Executes on button press in Stop_Recording.
function Stop_Recording_Callback(hObject, eventdata, handles)
% hObject    handle to Stop_Recording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);
myGUIdata.audioData = getaudiodata(myGUIdata.recordObj1);
%disp('timer ends')
%set(myGUIdata.startButton,'enable','off');

%off button
set(myGUIdata.Stop_Recording, 'enable', 'off');
%on button
set(myGUIdata.Save_RecordingSound, 'enable', 'on');
set(myGUIdata.Start_Recording, 'enable', 'on');
set(myGUIdata.Play_Original_Sound, 'enable', 'on');
set(myGUIdata.make_Improve_sound, 'enable', 'on');


%{
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms)
    case 'off'
    otherwise
        disp('timer is bloken!');
end;
%}

stop(myGUIdata.recordObj1);
guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end

% --- Executes on button press in Save_RecordingSound.
function Save_RecordingSound_Callback(hObject, eventdata, handles)
% hObject    handle to Save_RecordingSound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

set(myGUIdata.Save_RecordingSound, 'enable', 'off');

outFileName = ['audioIn' datestr(now,30) '.wav'];
[file,path] = uiputfile(outFileName,'Save the captured data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;
%wavwrite(myGUIdata.audioData,myGUIdata.samplingFrequency,16,[path file]);
audiowrite([path file],myGUIdata.audioData,myGUIdata.samplingFrequency,'BitsPerSample',16);
end

% --- Executes on button press in make_Improve_sound.
function make_Improve_sound_Callback(hObject, eventdata, handles)
% hObject    handle to make_Improve_sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

%button processSave_RecordingSound, 'enable', 'off');
set(myGUIdata.make_Improve_sound, 'enable', 'off');

set(myGUIdata.Save_RecordingSound, 'enable', 'off');
set(myGUIdata.Save_Modified_Sound, 'enable', 'off');
set(myGUIdata.Play_modified_Sound, 'enable', 'off');
set(myGUIdata.Play_Original_Sound, 'enable', 'off');

pause(1)

x = myGUIdata.audioData;
fs = myGUIdata.samplingFrequency;

crSndStr = TestgenerateIndividualDictionaryImproveSound(x,fs);
myGUIdata.crSndStr = crSndStr;

dictionaryDirPath = 'manipulationDictionary2014122401';
dictionaryDir = dir(dictionaryDirPath);
myGUIdata.dictionaryDir = dictionaryDir;

cmpVar = 0;

for dInd = 1:length(dictionaryDir)
    if strfind(dictionaryDir(dInd).name,'.mat') > 0
        dictionary = load([dictionaryDirPath '/' dictionaryDir(dInd).name]);        
        mStr = dictionary.manipulationStructure;
        crDiffDis = (crSndStr.tmpPrm.tmpDistance - mStr.tmpPrm.DistUnAtt)^2;
        if crDiffDis > cmpVar
            cmpVar = crDiffDis;
            minInd = dInd;
            myGUIdata.minInd = minInd;
        end
    end
end

% improve voice processing
dictionary = load([dictionaryDirPath '/' dictionaryDir(minInd).name]);
crImpSoundStr = generateIndividualDictionaryImproveSound(x,fs,dictionary);

myGUIdata.crImpSoundStr = crImpSoundStr;

%on button
set(myGUIdata.Save_RecordingSound, 'enable', 'on');
set(myGUIdata.Save_Modified_Sound, 'enable', 'on');
set(myGUIdata.Start_Recording, 'enable', 'on');
set(myGUIdata.Play_modified_Sound, 'enable', 'on');
set(myGUIdata.Play_Original_Sound, 'enable', 'on');

guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end

% --- Executes on button press in Play_modified_Sound.
function Play_modified_Sound_Callback(hObject, eventdata, handles)
% hObject    handle to Play_modified_Sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

modSoundStr = myGUIdata.crImpSoundStr;

x = modSoundStr.all.outputBuffer;
fs = modSoundStr.lpcStruct.Fs;
minus4 = max(1,round((length(x)/fs-4)*fs));
sound(x(minus4:end)/max(abs(x(minus4:end)))*0.99,fs);

end

% --- Executes on button press in Play_Original_Sound.
function Play_Original_Sound_Callback(hObject, eventdata, handles)
% hObject    handle to Play_Original_Sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

x = myGUIdata.audioData;
fs = myGUIdata.samplingFrequency;
minus4 = max(1,round((length(x)/fs-4)*fs));
sound(x(minus4:end)/max(abs(x(minus4:end)))*0.99,fs);

end

% --- Executes on button press in Save_Modified_Sound.
function Save_Modified_Sound_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Modified_Sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

set(myGUIdata.Save_Modified_Sound, 'enable', 'off');

modSoundStr = myGUIdata.crImpSoundStr;

x = modSoundStr.all.outputBuffer;
fs = modSoundStr.lpcStruct.Fs;

outFileName = ['modifiedAudioIn' datestr(now,30) '.wav'];
[file,path] = uiputfile(outFileName,'Save the captured data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;
%wavwrite(myGUIdata.audioData,myGUIdata.samplingFrequency,16,[path file]);
audiowrite([path file],x,fs,'BitsPerSample',16);

end


% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

%disp('timer ends')
%{
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms)
end;
%}
stop(myGUIdata.recordObj1);

delete(myGUIdata.recordObj1);

close(handles.AttractivenessImprovementMainGUI);
end
