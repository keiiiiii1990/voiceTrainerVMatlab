function varargout = vtlEstimator(varargin)
%   Realtime FFT analyzer with waveform sub-display. Type:
%   vtlEstimator
%   to start.

%   Designed and coded by Hideki Kawahara (kawahara AT sys.wakayama-u.ac.jp)
%   28/Nov./2013
%   29/Nov./2013 minor bug fix
%   No warranty:
%   This is a sample code. Use this as you like.

% VTLESTIMATOR MATLAB code for vtlEstimator.fig
%      VTLESTIMATOR, by itself, creates a new VTLESTIMATOR or raises the existing
%      singleton*.
%
%      H = VTLESTIMATOR returns the handle to a new VTLESTIMATOR or the handle to
%      the existing singleton*.
%
%      VTLESTIMATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VTLESTIMATOR.M with the given input arguments.
%
%      VTLESTIMATOR('Property','Value',...) creates a new VTLESTIMATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vtlEstimator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vtlEstimator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vtlEstimator

% Last Modified by GUIDE v2.5 12-Dec-2013 15:55:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @vtlEstimator_OpeningFcn, ...
    'gui_OutputFcn',  @vtlEstimator_OutputFcn, ...
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


% --- Executes just before vtlEstimator is made visible.
function vtlEstimator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vtlEstimator (see VARARGIN)

global handleForTimer;
% Choose default command line output for vtlEstimator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%handles
clear global vtlEstimatorGlobal
delete(timerfindall);
%pause(0.05)
initializeDisplay(handles);
myGUIdata = guidata(handles.multiScopeMainGUI);
%set(myGUIdata.largeViewTypePopup,'enable','off');
set(myGUIdata.vtlDisplay,'visible','off');
set(myGUIdata.saveButton,'enable','off');
set(myGUIdata.vtlEstimateButton,'enable','off');
%pause(0.5)
timerEventInterval = 0.1; % in second
timer50ms = timer('TimerFcn',@synchDrawGUI, 'Period',timerEventInterval,'ExecutionMode','fixedRate');
handleForTimer = handles.multiScopeMainGUI; % global
myGUIdata.timer50ms = timer50ms;
myGUIdata.smallviewerWidth = 30; % 30 ms is defaule
myGUIdata.samplingFrequency = 44100;
myGUIdata.recordObj1 = audiorecorder(myGUIdata.samplingFrequency,24,1);
record(myGUIdata.recordObj1)
switch get(myGUIdata.recordObj1,'Running')
    case 'on'
        stop(myGUIdata.recordObj1);
end
guidata(handles.multiScopeMainGUI,myGUIdata);
startButton_Callback(hObject, eventdata, handles)
% UIWAIT makes vtlEstimator wait for user response (see UIRESUME)
% uiwait(handles.multiScopeMainGUI);
end

function initializeDisplay(handles)
myGUIdata = guidata(handles.multiScopeMainGUI);
myGUIdata.maxAudioRecorderCount = 200;
myGUIdata.audioRecorderCount = myGUIdata.maxAudioRecorderCount;
myGUIdata.maxLevelIndicator = -100*ones(myGUIdata.maxAudioRecorderCount,1);
myGUIdata.yMax = 1;
axes(myGUIdata.smallViewerAxis);
myGUIdata.smallViewerPlotHandle = plot(randn(1000,1));
set(myGUIdata.smallViewerAxis,'xtick',[],'ytick',[]);
axes(myGUIdata.largeViewerAxis);
fs = 44100;
dataLength = round(30/1000*fs);
fftl = 2.0.^ceil(log2(dataLength));
fAxis = (0:fftl-1)/fftl*fs;
w = blackman(dataLength);
pw = 20*log10(abs(fft(randn(dataLength,1).*w,fftl)/sqrt(sum(w.^2))));
myGUIdata.axisType = 'Logarithmic';
myGUIdata.window = w;
myGUIdata.fAxis = fAxis;
switch myGUIdata.axisType
    case 'Linear'
        myGUIdata.largeViewerPlotHandle = plot(fAxis,pw);grid on;
        axis([0 fs/2 [-90 20]]);
    case 'Logarithmic'
        myGUIdata.largeViewerPlotHandle = semilogx(fAxis,pw);grid on;
        axis([10 fs/2 [-90 20]]);
end;
set(gca,'fontsize',15);
xlabel('frequency (Hz)');
ylabel('level (dB)');
axes(myGUIdata.wholeViewerAxis);
myGUIdata.wholeViewerHandle = plot(myGUIdata.maxLevelIndicator);
axis([0 myGUIdata.maxAudioRecorderCount -100 0]);
set(myGUIdata.wholeViewerAxis,'xtick',[],'ylim',[-80 0]);grid on;
guidata(handles.multiScopeMainGUI,myGUIdata);
end

function synchDrawGUI(obj, event, string_arg)
global handleForTimer;
myGUIdata = guidata(handleForTimer);
myGUIdata.smallviewerWidth = get(myGUIdata.radiobutton10ms,'value')*10+ ...
    get(myGUIdata.radiobutton30ms,'value')*30+ ...
    get(myGUIdata.radiobutton100ms,'value')*100+ ...
    get(myGUIdata.radiobutton300ms,'value')*300;
numberOfSamples = round(myGUIdata.smallviewerWidth*myGUIdata.samplingFrequency/1000);
if get(myGUIdata.recordObj1,'TotalSamples') > numberOfSamples
tmpAudio = getaudiodata(myGUIdata.recordObj1);
currentPoint = length(tmpAudio);
xdata = 1:numberOfSamples;
fs = myGUIdata.samplingFrequency;
%disp(myGUIdata.audioRecorderCount)
if length(currentPoint-numberOfSamples+1:currentPoint) > 10
    ydata = tmpAudio(currentPoint-numberOfSamples+1:currentPoint);
    myGUIdata.audioRecorderCount = myGUIdata.audioRecorderCount-1;
    set(myGUIdata.smallViewerPlotHandle,'xdata',xdata,'ydata',ydata);
    if myGUIdata.yMax < max(abs(ydata))
        myGUIdata.yMax = max(abs(ydata));
    else
        myGUIdata.yMax = myGUIdata.yMax*0.8;
    end;
    set(myGUIdata.smallViewerAxis,'xlim',[0 numberOfSamples],'ylim',myGUIdata.yMax*[-1 1]);
    fftl = 2^ceil(log2(numberOfSamples));
    fAxis = (0:fftl-1)/fftl*fs;
    switch get(myGUIdata.windowTypePopup,'value')
        case 1
            w = blackman(numberOfSamples);
        case 2
            w = hamming(numberOfSamples);
        case 3
            w = hanning(numberOfSamples);
        case 4
            w = bartlett(numberOfSamples);
        case 5
            w = ones(numberOfSamples,1);
        case 6
            w = nuttallwin(numberOfSamples);
        otherwise
            w = blackman(numberOfSamples);
    end;
    pw = 20*log10(abs(fft(ydata.*w,fftl)/sqrt(sum(w.^2))));
    set(myGUIdata.largeViewerPlotHandle,'xdata',fAxis,'ydata',pw);
    myGUIdata.maxLevelIndicator(max(1,myGUIdata.maxAudioRecorderCount-myGUIdata.audioRecorderCount)) ...
        = max(20*log10(abs(ydata)));
    set(myGUIdata.wholeViewerHandle,'ydata',myGUIdata.maxLevelIndicator);
else
    disp('overrun!')
end;
if myGUIdata.audioRecorderCount < 0
    switch get(myGUIdata.timer50ms,'running')
        case 'on'
            stop(myGUIdata.timer50ms);
    end
    stop(myGUIdata.recordObj1);
    record(myGUIdata.recordObj1);
    myGUIdata.audioRecorderCount = myGUIdata.maxAudioRecorderCount;
    myGUIdata.maxLevelIndicator = 0;
    %    switch get(myGUIdata.timer50ms,'running')
    %        case 'on'
    %            stop(myGUIdata.timer50ms)
    %    end;
    switch get(myGUIdata.timer50ms,'running')
        case 'off'
            start(myGUIdata.timer50ms);
    end
end;
guidata(handleForTimer,myGUIdata);
else
    disp(['Recorded data is not enough! Skipping this interruption....at ' datestr(now,30)]);
end;
end

% --- Outputs from this function are returned to the command line.
function varargout = vtlEstimator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global handleForTimer
myGUIdata = guidata(handles.multiScopeMainGUI);
set(myGUIdata.saveButton,'enable','off');
set(myGUIdata.startButton,'enable','off');
set(myGUIdata.stopButton,'enable','on');
set(myGUIdata.vtlDisplay,'visible','off');
set(myGUIdata.vtlEstimateButton,'enable','off');
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms);
end
myGUIdata.audioRecorderCount = myGUIdata.maxAudioRecorderCount;
myGUIdata.maxLevelIndicator = -100*ones(myGUIdata.maxAudioRecorderCount,1);
myGUIdata.yMax = 1;
record(myGUIdata.recordObj1);
switch get(myGUIdata.timer50ms,'running')
    case 'off'
        start(myGUIdata.timer50ms);
    case 'on'
    otherwise
        disp('timer is bloken!');
end
guidata(handles.multiScopeMainGUI,myGUIdata);
end


% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, eventdata, handles)
% hObject    handle to stopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.multiScopeMainGUI);
myGUIdata.audioData = getaudiodata(myGUIdata.recordObj1);
%disp('timer ends')
%set(myGUIdata.startButton,'enable','off');
set(myGUIdata.saveButton,'enable','on');
set(myGUIdata.startButton,'enable','on');
set(myGUIdata.stopButton,'enable','off');
set(myGUIdata.vtlEstimateButton,'enable','on');
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms)
    case 'off'
    otherwise
        disp('timer is bloken!');
end;
stop(myGUIdata.recordObj1);
guidata(handles.multiScopeMainGUI,myGUIdata);
end


% --- Executes on selection change in windowTypePopup.
function windowTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to windowTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns windowTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from windowTypePopup
end


% --- Executes during object creation, after setting all properties.
function windowTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in radiobutton10ms.
function radiobutton10ms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10ms
myGUIdata = guidata(handles.multiScopeMainGUI);
set(handles.radiobutton10ms,'value',1);
set(handles.radiobutton30ms,'value',0);
set(handles.radiobutton100ms,'value',0);
set(handles.radiobutton300ms,'value',0);
end


% --- Executes on button press in radiobutton30ms.
function radiobutton30ms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton30ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton30ms
myGUIdata = guidata(handles.multiScopeMainGUI);
set(handles.radiobutton10ms,'value',0);
set(handles.radiobutton30ms,'value',1);
set(handles.radiobutton100ms,'value',0);
set(handles.radiobutton300ms,'value',0);
end


% --- Executes on button press in radiobutton100ms.
function radiobutton100ms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton100ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton100ms
myGUIdata = guidata(handles.multiScopeMainGUI);
set(handles.radiobutton10ms,'value',0);
set(handles.radiobutton30ms,'value',0);
set(handles.radiobutton100ms,'value',1);
set(handles.radiobutton300ms,'value',0);
end


% --- Executes on button press in radiobutton300ms.
function radiobutton300ms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton300ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton300ms
myGUIdata = guidata(handles.multiScopeMainGUI);
set(handles.radiobutton10ms,'value',0);
set(handles.radiobutton30ms,'value',0);
set(handles.radiobutton100ms,'value',0);
set(handles.radiobutton300ms,'value',1);
end


% --- Executes on selection change in axisTypePopup.
function axisTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to axisTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns axisTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from axisTypePopup
myGUIdata = guidata(handles.multiScopeMainGUI);
fs = myGUIdata.samplingFrequency;
switch get(myGUIdata.axisTypePopup,'value')
    case 1
        set(myGUIdata.largeViewerAxis,'xlim',[10 fs/2],'xscale','log');
    case 2
        set(myGUIdata.largeViewerAxis,'xlim',[0 fs/2],'xscale','linear');
    case 3
        set(myGUIdata.largeViewerAxis,'xlim',[0 10000],'xscale','linear');
    case 4
        set(myGUIdata.largeViewerAxis,'xlim',[0 8000],'xscale','linear');
    case 5
        set(myGUIdata.largeViewerAxis,'xlim',[0 4000],'xscale','linear');
    otherwise
        set(myGUIdata.largeViewerAxis,'xlim',[10 fs/2],'xscale','log');
end;
end


% --- Executes during object creation, after setting all properties.
function axisTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axisTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.multiScopeMainGUI);
%disp('timer ends')
switch get(myGUIdata.timer50ms,'running')
    case 'on'
        stop(myGUIdata.timer50ms)
end;
stop(myGUIdata.recordObj1);
delete(myGUIdata.timer50ms);
delete(myGUIdata.recordObj1);
close(handles.multiScopeMainGUI);
end


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.multiScopeMainGUI);
outFileName = ['audioIn' datestr(now,30) '.wav'];
[file,path] = uiputfile(outFileName,'Save the captured data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;
wavwrite(myGUIdata.audioData,myGUIdata.samplingFrequency,16,[path file]);
end


% --- Executes on button press in vtlEstimateButton.
function vtlEstimateButton_Callback(hObject, eventdata, handles)
% hObject    handle to vtlEstimateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.multiScopeMainGUI);
fs = myGUIdata.samplingFrequency;
set(handles.multiScopeMainGUI,'pointer','watch');drawnow
x = hanningHPF(myGUIdata.audioData,fs,70);
vtlStructure = textIndependentVTLKalman(x,fs);
set(myGUIdata.vtlDisplay,'string',[num2str(vtlStructure.vtlBestInCm,'%3.2f') ' cm']);
set(myGUIdata.vtlDisplay,'visible','on');
set(handles.multiScopeMainGUI,'pointer','arrow');drawnow

end