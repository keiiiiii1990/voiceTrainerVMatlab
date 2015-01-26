function varargout = attractiveImproverV03(varargin)
% ATTRACTIVEIMPROVERV03 MATLAB code for attractiveImproverV03.fig
%      ATTRACTIVEIMPROVERV03, by itself, creates a new ATTRACTIVEIMPROVERV03 or raises the existing
%      singleton*.
%
%      H = ATTRACTIVEIMPROVERV03 returns the handle to a new ATTRACTIVEIMPROVERV03 or the handle to
%      the existing singleton*.
%
%      ATTRACTIVEIMPROVERV03('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ATTRACTIVEIMPROVERV03.M with the given input arguments.
%
%      ATTRACTIVEIMPROVERV03('Property','Value',...) creates a new ATTRACTIVEIMPROVERV03 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before attractiveImproverV03_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to attractiveImproverV03_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help attractiveImproverV03

% Last Modified by GUIDE v2.5 06-Jan-2015 10:40:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @attractiveImproverV03_OpeningFcn, ...
                   'gui_OutputFcn',  @attractiveImproverV03_OutputFcn, ...
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

% --- Executes just before attractiveImproverV03 is made visible.
function attractiveImproverV03_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to attractiveImproverV03 (see VARARGIN)

% Choose default command line output for attractiveImproverV03
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%% mycode
initializeDisplay(handles);

% UIWAIT makes attractiveImproverV03 wait for user response (see UIRESUME)
% uiwait(handles.AttractivenessImprovementMainGUI);

%%%%%%%%%%% my code
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

%For Initialize set parameter
set(myGUIdata.OrgVTAFDispBtn,'value',0);
set(myGUIdata.ModVTAFDispBtn,'value',0);
initializeDisplay(handles);

set(myGUIdata.Start_Recording, 'enable', 'on');

set(myGUIdata.Stop_Recording, 'enable', 'off');
set(myGUIdata.Save_RecordingSound, 'enable', 'off');
set(myGUIdata.Play_modified_Sound, 'enable', 'off');
set(myGUIdata.Play_Original_Sound, 'enable', 'off');
set(myGUIdata.Save_Modified_Sound, 'enable', 'off');
set(myGUIdata.make_Improve_sound, 'enable', 'off');
set(myGUIdata.StartMov_OrgVTAF, 'enable', 'off');
set(myGUIdata.StartMov_ModVTAF, 'enable', 'off');

set(myGUIdata.OrgVTAFDispBtn,'enable','off');
set(myGUIdata.ModVTAFDispBtn,'enable','off');

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
function varargout = attractiveImproverV03_OutputFcn(hObject, eventdata, handles) 
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
myGUIdata.player = audioplayer(randn(1000,1),fs);
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
myGUIdata.tractRotateOrgVTAF = rotate3d;
set(myGUIdata.tractRotateOrgVTAF,'enable','on');

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
%rotate3d on;
myGUIdata.tractRotateModVTAF = rotate3d;
set(myGUIdata.tractRotateModVTAF,'enable','on');

if get(myGUIdata.OrgVTAFDispBtn,'value') == 1
    set(myGUIdata.tract3DOrg,'visible','on');
    set(myGUIdata.tractRotateOrgVTAF,'enable','on');
    
elseif get(myGUIdata.OrgVTAFDispBtn,'value') == 0
    set(myGUIdata.tract3DOrg,'visible','off');
    set(myGUIdata.tractRotateOrgVTAF,'enable','off');
end

if get(myGUIdata.ModVTAFDispBtn,'value') == 1
    set(myGUIdata.tract3DMod,'visible','on');
    set(myGUIdata.tractRotateModVTAF,'enable','on');
elseif get(myGUIdata.ModVTAFDispBtn,'value') == 0
    set(myGUIdata.tract3DMod,'visible','off');
    set(myGUIdata.tractRotateModVTAF,'enable','off');
end

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

set(myGUIdata.OrgVTAFDispBtn,'enable','off');
set(myGUIdata.OrgVTAFDispBtn,'value',0);
OrgVTAFDispBtn_Callback(hObject, eventdata, handles);
set(myGUIdata.ModVTAFDispBtn,'enable','off');
set(myGUIdata.ModVTAFDispBtn,'value',0);
ModVTAFDispBtn_Callback(hObject, eventdata, handles);

%off uiController
set(handles.allVTAFFrmSlider,'enable','off');
set(handles.OrgVTAFFrmSlider,'enable','off');
set(handles.ModVTAFFrmSlider,'enable','off');

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
set(myGUIdata.OrgVTAFDispBtn,'enable','on');
set(myGUIdata.OrgVTAFDispBtn,'value',1);
OrgVTAFDispBtn_Callback(hObject, eventdata, handles);

set(myGUIdata.Save_RecordingSound, 'enable', 'on');
set(myGUIdata.Start_Recording, 'enable', 'on');
set(myGUIdata.Play_Original_Sound, 'enable', 'on');
set(myGUIdata.StartMov_OrgVTAF, 'enable', 'on');
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


%calculate vtaf
x = myGUIdata.audioData;
fs = myGUIdata.samplingFrequency;

%%f0 relation paramter
f0Struct = exF0candidatesTSTRAIGHTGB(x,fs);
voicedF0 = f0Struct.f0(f0Struct.periodicityLevel>0.5);
unVoicedF0 = f0Struct.f0(f0Struct.periodicityLevel<0.5);

ax = x;samplingFrequency = fs;
VtafStruct = struct;
VtafStruct.ax = ax;
VtafStruct.samplingFrequency = samplingFrequency;

frameLengthInMs = 25;
frameShiftInMs = 5;

responseLengthInMs = 20;
Fs = 8000;
convertedStr = convertTo8kHzSampledSignal(ax,samplingFrequency);
lpcStructure = plainLPCAnalysis(convertedStr.signal,convertedStr.samplingFrequency,frameShiftInMs);
VtafStruct.lpcStructure = lpcStructure;
f0Struct = higherSymKalmanWithTIFupdate(ax,samplingFrequency);
medianF0 = median(f0Struct.F0(f0Struct.latentSDcoeff<1.07));
signalTime = (0:1/Fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);

for ii = 1:Fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/Fs*2*pi)).*(f0InSignalTime*ii<Fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/Fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*Fs):round(lpcStructure.frameShiftInMs/1000*Fs))';
responseLength = round(responseLengthInMs/1000*Fs);
fftLSynth = 2^ceil(log2(responseLength+length(baseIndex)+1));
olaIndex = (1:fftLSynth)';
for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStructure.preEmphasis*cos(theta(:))) ...
        +lpcStructure.cepstrumList(ii,1)*cos(theta(:)) ...
        +lpcStructure.cepstrumList(ii,2)*cos(2*theta(:));
    LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStructure.powerList(ii)/sqrt(f0Struct.F0(ii));
    cepstrum = ifft(log(fixedLPCspectrumInPower)/2);
    cepstrum(fftl/2:fftl) = 0;
    cepstrum(2:fftl/2) = cepstrum(2:fftl/2)*2;
    minimumPhaseSpectrum = exp(fft(cepstrum));
    minimumPhaseResponse = real(ifft(minimumPhaseSpectrum));
    segment = excitation(max(1,min(length(excitation),baseIndex+round(lpcStructure.temporalPosition(ii)*8000))));
    filterResponse = minimumPhaseResponse(1:responseLength);
    response = real(ifft(fft(filterResponse,fftLSynth).*fft(segment,fftLSynth)));
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*Fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*Fs)))+response;
end;
waveform = outputBuffer;
Fs = 8000;
VtafStruct.waveform = waveform;
VtafStruct.Fs = Fs;
myGUIdata.VtafStruct = VtafStruct;


set(handles.OrgVTAFFrmSlider,'enable','on');
[mv,minIndex] = min(VtafStruct.lpcStructure.temporalPosition(floor(length(VtafStruct.lpcStructure.temporalPosition)/2):end));

Max = (length(lpcStructure.temporalPosition(:,1)));
majorS = ceil(length(lpcStructure.temporalPosition(:,1))/10);
Value = minIndex;

%------ information for waveform monitor (lissajous handle)
lowpassLength = round(1.7*myGUIdata.samplingFrequency/medianF0/2)*2+1;
blackmanLowPass = blackman(lowpassLength);
xxLPF = fftfilt(diff(blackmanLowPass/sum(blackmanLowPass))*100,[myGUIdata.audioData;zeros(lowpassLength,1)]);
myGUIdata.fundamentalComponent = xxLPF((1:length(myGUIdata.audioData))+floor(lowpassLength/2));

set(handles.OrgVTAFFrmSlider,'Max',Max,'SliderStep',[1/majorS 1/10],'Value',Value);

guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end

%%private function-- crossSection = tubeDisplay(logArea)
%this use function is...
% stop_recording_callback()
function crossSection = tubeDisplay(logArea)
areaList = exp(logArea);
crossSection = sqrt(areaList/sum(areaList));
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
crImpSoundStr = generateIndividualDictionaryImproveSoundV02(x,fs,dictionary);

myGUIdata.crImpSoundStr = crImpSoundStr;

%on button
set(myGUIdata.Save_RecordingSound, 'enable', 'on');
set(myGUIdata.Save_Modified_Sound, 'enable', 'on');
set(myGUIdata.Start_Recording, 'enable', 'on');
set(myGUIdata.Play_modified_Sound, 'enable', 'on');
set(myGUIdata.Play_Original_Sound, 'enable', 'on');
set(myGUIdata.StartMov_ModVTAF, 'enable', 'on');

set(myGUIdata.tractRotateModVTAF,'enable','on');
set(myGUIdata.ModVTAFDispBtn,'enable','on');
set(myGUIdata.ModVTAFDispBtn,'value',1);
ModVTAFDispBtn_Callback(hObject, eventdata, handles);

set(handles.allVTAFFrmSlider,'enable','on');
set(handles.ModVTAFFrmSlider,'enable','on');

[mv,minIndex] = min(crImpSoundStr.newLpcStructure.lpcStructure.temporalPosition(floor(length(crImpSoundStr.newLpcStructure.lpcStructure.temporalPosition)/2):end));
Max = (length(crImpSoundStr.newLpcStructure.lpcStructure.temporalPosition(:,1)));
majorS = ceil(length(crImpSoundStr.newLpcStructure.lpcStructure.temporalPosition(:,1))/10);
%shO.SliderStep = [1/majorS 1/10];
Value = minIndex;

set(handles.ModVTAFFrmSlider,'Max',Max,'SliderStep',[1/majorS 1/10],'Value',Value);
set(handles.allVTAFFrmSlider,'Max',Max,'SliderStep',[1/majorS 1/10],'Value',Value);

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
VtafStruct = myGUIdata.VtafStruct ;

minus4 = max(1,round((length(x)/fs-4)*fs));
%sound(x(minus4:end)/max(abs(x(minus4:end)))*0.99,fs);
%{
myGUIdata.minus4 = minus4;
myGUIdata.audoPlayTimerPeriod = 0.057;
switch get(myGUIdata.player,'running')
    case 'off'
        myGUIdata.player = audioplayer(x(minus4:end)/max(abs(x(minus4:end)))*0.9,fs);
        set(myGUIdata.Play_Original_Sound,'string','STOP');
        guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
        %set(myGUIdata.sgramCursorHandle,'visible','on');
        %******************************           Importance part          *******************************%
        set(myGUIdata.player,'userdata',handles.AttractivenessImprovementMainGUI,'TimerFcn', ...
            @sliderCursorUpdateTimerFunction,'TimerPeriod',myGUIdata.audoPlayTimerPeriod);
        
        %set(myGUIdata.counterText,'visible','on');
        %set(myGUIdata.lissajousPlotHandle,'visible','on');
        playblocking(myGUIdata.player);
        set(myGUIdata.Play_Original_Sound,'string','PLAY');
        %set(myGUIdata.counterText,'visible','off');
        %set(myGUIdata.lissajousPlotHandle,'visible','off');
        %set(myGUIdata.sgramCursorHandle,'visible','off');
    case 'on'
        stop(myGUIdata.player);
        set(myGUIdata.Play_Original_Sound,'string','PLAY');
        %set(myGUIdata.counterText,'visible','off');
        %set(myGUIdata.lissajousPlotHandle,'visible','off');
        %set(myGUIdata.sgramCursorHandle,'visible','off');
end;
%}

%
%axes(handles.orgVTAFView);
%vLen = length(VtafStruct.lpcStructure.logAreaMatrix(:,1));
%Ssec = VtafStruct.lpcStructure.temporalPosition(end);
%fpSec = 
%set(myGUIdata.textForDebug01,'String',num2str(vLen));
%{
for ii = 1:vLen
    logArea = VtafStruct.lpcStructure.logAreaMatrix(ii,:);
    cSec = tubeDisplay(logArea(:));
    [X,Y,Z] = cylinder(cSec,40);
    set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
    drawnow
    %M(ii) = getframe(gca);
end
%}
%movie(M,1,30);

%{
for ii = 1:vLen
    logArea = VtafStruct.lpcStructure.logAreaMatrix(ii,:);
    cSec = tubeDisplay(logArea(:));
    [X,Y,Z] = cylinder(cSec,40);
    set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
    drawnow
end
%}

end

%praivate function
function sliderCursorUpdateTimerFunction(obj, event, string_arg)
%obj is "audioplayer"Object , which has "userData" property.

handleForTimer = get(obj,'userData');
myGUIdata = guidata(handleForTimer);

VtafStruct = myGUIdata.VtafStruct ;

currentPoint = get(myGUIdata.player,'CurrentSample');
fs = get(myGUIdata.player,'SampleRate');
index30ms = 1:round(0.03*fs);
TotalSamples = get(myGUIdata.player,'TotalSamples');

xx = myGUIdata.fundamentalComponent(min(TotalSamples,currentPoint+index30ms));

triggerPoint = find(xx(1:end-1).*xx(2:end)<0 & xx(1:end-1)<0, 1 );
ydata = myGUIdata.audioData(min(TotalSamples,currentPoint+index30ms+triggerPoint));

axes(handles.orgVTAFView);

logArea = VtafStruct.lpcStructure.logAreaMatrix(ii,:);
cSec = tubeDisplay(logArea(:));
[X,Y,Z] = cylinder(cSec,40);

if sum(isnan(ydata))>0 || length(ydata) <1 || max(abs(ydata)) == 0
%    set(myGUIdata.lissajousAxis,'ylim',[-1 1],'xlim',[1 index30ms(end)]);
    set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
else
%    set(myGUIdata.lissajousAxis,'ylim',max(abs(ydata))*[-1 1],'xlim',[1 index30ms(end)]);
end;

%set(myGUIdata.sgramCursorHandle,'xdata',currentPoint/fs*[1 1]);

end

%when record sound play, update monitor function
function monitorUpdateTimerFunction(obj, event, string_arg)
%   timer function for updating display

%   14/May/2014
%   by Hideki Kawahara
handleForTimer = get(obj,'userData');
%myGUIdata = get(handleForTimer,'userdata');
myGUIdata = guidata(handleForTimer);
currentPoint = get(myGUIdata.recordObj1,'CurrentSample');
fs = get(myGUIdata.recordObj1,'SampleRate');
index30ms = -round(0.03*fs):0;
%TotalSamples = get(myGUIdata.recorder,'TotalSamples');
y = getaudiodata(myGUIdata.recorder);
ydata = y(max(1,min(length(y),length(y)+index30ms)));
set(myGUIdata.lissajousAxis,'ylim',max(abs(ydata))*[-1 1],'xlim',[index30ms(1) index30ms(end)]);
set(myGUIdata.lissajousPlotHandle,'xdata',index30ms,'ydata',ydata);
set(myGUIdata.counterText,'string',[num2str(currentPoint/fs*1000,'%05.0f') '  ms']);
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


% --- Executes on button press in OrgVTAFDispBtn.
function OrgVTAFDispBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OrgVTAFDispBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OrgVTAFDispBtn
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

axes(handles.orgVTAFView);
if get(myGUIdata.OrgVTAFDispBtn,'value') == 1
    %set(myGUIdata.tract3DOrg,'visible','on')
    set(get(gca,'children'),'visible','on');
    set(myGUIdata.tractRotateOrgVTAF,'enable','on');
elseif get(myGUIdata.OrgVTAFDispBtn,'value') == 0
    %set(myGUIdata.tract3DOrg,'visible','off')
    set(get(gca,'children'),'visible','off');
    set(myGUIdata.tractRotateOrgVTAF,'enable','off');
end

guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end


% --- Executes on button press in ModVTAFDispBtn.
function ModVTAFDispBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ModVTAFDispBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ModVTAFDispBtn
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

axes(handles.modVTAFView);
if get(myGUIdata.ModVTAFDispBtn,'value') == 1
    %set(myGUIdata.tract3DOrg,'visible','on')
    set(get(gca,'children'),'visible','on');
    set(myGUIdata.tractRotateModVTAF,'enable','on');
elseif get(myGUIdata.ModVTAFDispBtn,'value') == 0
    %set(myGUIdata.tract3DOrg,'visible','off')
    set(get(gca,'children'),'visible','off');
    set(myGUIdata.tractRotateModVTAF,'enable','off');
end

guidata(handles.AttractivenessImprovementMainGUI,myGUIdata);
end


% --- Executes on slider movement.
function OrgVTAFFrmSlider_Callback(hObject, eventdata, handles)
% hObject    handle to OrgVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);
VtafStruct = myGUIdata.VtafStruct ;

axes(handles.orgVTAFView);
if floor(get(hObject,'Value')) > 0
VTAFFrmIndex = floor(get(hObject,'Value'));
logArea = VtafStruct.lpcStructure.logAreaMatrix(VTAFFrmIndex,:);
cSec = tubeDisplay(logArea(:));
[X,Y,Z] = cylinder(cSec,40);
set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
drawnow
end

end

% --- Executes during object creation, after setting all properties.
function OrgVTAFFrmSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OrgVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on slider movement.
function ModVTAFFrmSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ModVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

crImpSoundStr = myGUIdata.crImpSoundStr ;

axes(handles.modVTAFView);

if floor(get(hObject,'Value')) > 0
VTAFFrmIndex = floor(get(hObject,'Value'));
logArea = crImpSoundStr.newLpcStructure.newLogAreaMatrix(VTAFFrmIndex,:);
cSec = tubeDisplay(logArea(:));
[X,Y,Z] = cylinder(cSec,40);
set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
drawnow
end

end

% --- Executes during object creation, after setting all properties.
function ModVTAFFrmSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on slider movement.
function allVTAFFrmSlider_Callback(hObject, eventdata, handles)
% hObject    handle to allVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);

VtafStruct = myGUIdata.VtafStruct ;
crImpSoundStr = myGUIdata.crImpSoundStr ;

axes(handles.orgVTAFView);
oA = gca;
axes(handles.modVTAFView);
mA = gca;


if get(hObject,'Min') < get(hObject,'Value')&&  get(hObject,'Value') < get(hObject,'Max')
    if isequal(get(get(oA,'children'),'visible'),'on') && isequal(get(get(mA,'children'),'visible'),'on')
        VTAFFrmIndex = floor(get(hObject,'Value'));
        %org vtaf processing
        ologArea = VtafStruct.lpcStructure.logAreaMatrix(VTAFFrmIndex,:);
        ocSec = tubeDisplay(ologArea(:));
        [oX,oY,oZ] = cylinder(ocSec,40);
        set(get(oA,'children'),'xdata',oZ,'ydata',oY,'zdata',oX,'visible','on');
        %drawnow        
        %mod vtaf processing
        mlogArea = crImpSoundStr.newLpcStructure.newLogAreaMatrix(VTAFFrmIndex,:);
        mcSec = tubeDisplay(mlogArea(:));
        [mX,mY,mZ] = cylinder(mcSec,40);
        set(get(mA,'children'),'xdata',mZ,'ydata',mY,'zdata',mX,'visible','on');
        drawnow
        
        set(handles.OrgVTAFFrmSlider,'Value',get(hObject,'Value'));
        set(handles.ModVTAFFrmSlider,'Value',get(hObject,'Value'));
    end
end
end

% --- Executes during object creation, after setting all properties.
function allVTAFFrmSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allVTAFFrmSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on button press in StartMov_OrgVTAF.
function StartMov_OrgVTAF_Callback(hObject, eventdata, handles)
% hObject    handle to StartMov_OrgVTAF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);
VtafStruct = myGUIdata.VtafStruct ;

axes(handles.orgVTAFView);
vLen = length(VtafStruct.lpcStructure.logAreaMatrix(:,1));
Ssec = VtafStruct.lpcStructure.temporalPosition(end);
for ii = 1:vLen
    logArea = VtafStruct.lpcStructure.logAreaMatrix(ii,:);
    cSec = tubeDisplay(logArea(:));
    [X,Y,Z] = cylinder(cSec,40);
    set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
    drawnow
    %M(ii) = getframe(gca);
end

%movie(M,1,30);


end


% --- Executes on button press in StartMov_ModVTAF.
function StartMov_ModVTAF_Callback(hObject, eventdata, handles)
% hObject    handle to StartMov_ModVTAF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.AttractivenessImprovementMainGUI);
crImpSoundStr = myGUIdata.crImpSoundStr ;

axes(handles.modVTAFView);
vLen = length(crImpSoundStr.newLpcStructure.newLogAreaMatrix(:,1));
Ssec = crImpSoundStr.newLpcStructure.lpcStructure.temporalPosition(end);
for ii = 1:vLen
    logArea = crImpSoundStr.newLpcStructure.newLogAreaMatrix(ii,:);
    cSec = tubeDisplay(logArea(:));
    [X,Y,Z] = cylinder(cSec,40);
    set(get(gca,'children'),'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
    drawnow
    %M(ii) = getframe(gca);
end

end