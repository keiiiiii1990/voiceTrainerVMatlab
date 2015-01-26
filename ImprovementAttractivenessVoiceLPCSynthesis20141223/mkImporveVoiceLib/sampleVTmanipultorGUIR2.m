function varargout = sampleVTmanipultorGUIR2(varargin)
% SAMPLEVTMANIPULTORGUIR2 MATLAB code for sampleVTmanipultorGUIR2.fig
%      SAMPLEVTMANIPULTORGUIR2, by itself, creates a new SAMPLEVTMANIPULTORGUIR2 or raises the existing
%      singleton*.
%
%      H = SAMPLEVTMANIPULTORGUIR2 returns the handle to a new SAMPLEVTMANIPULTORGUIR2 or the handle to
%      the existing singleton*.
%
%      SAMPLEVTMANIPULTORGUIR2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAMPLEVTMANIPULTORGUIR2.M with the given input arguments.
%
%      SAMPLEVTMANIPULTORGUIR2('Property','Value',...) creates a new SAMPLEVTMANIPULTORGUIR2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sampleVTmanipultorGUIR2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sampleVTmanipultorGUIR2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sampleVTmanipultorGUIR2
%   Designed and coded by Hideki Kawahara
%   09/June/2014
%   10/June/2014 first preliminary release
%   12/June/2014 revised display
%   13/June/2014 parameter save is enabled
%   14/June/2014 F0 and VTL control is added
%   15/June/2014 minor tuning
%   18/June/2014 minor bug fix
%   23/June/2014 made Windows compatible (hopefully)
%   17/July/2014 direct modifier control is added

% Last Modified by GUIDE v2.5 29-Jul-2014 10:44:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sampleVTmanipultorGUIR2_OpeningFcn, ...
    'gui_OutputFcn',  @sampleVTmanipultorGUIR2_OutputFcn, ...
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

% --- Executes just before sampleVTmanipultorGUIR2 is made visible.
function sampleVTmanipultorGUIR2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sampleVTmanipultorGUIR2 (see VARARGIN)

% Choose default command line output for sampleVTmanipultorGUIR2
handles.output = hObject;
guidata(hObject, handles);
myGUIdata = guidata(handles.VTmainpulatorGUI);

myGUIdata = initializeGUI(myGUIdata);
% Update handles structure
guidata(hObject, myGUIdata);

% UIWAIT makes sampleVTmanipultorGUIR2 wait for user response (see UIRESUME)
% uiwait(handles.VTmainpulatorGUI);
end

% --- Outputs from this function are returned to the command line.
function varargout = sampleVTmanipultorGUIR2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%--- private
function myGUIdata = initializeGUI(myGUIdata)
myGUIdata.debugmode = 'off';
if isequal(myGUIdata.debugmode,'on')
    clc
end
fs = 44100;
myGUIdata.samplingFrequency = fs;
axes(myGUIdata.lissajousAxis);
myGUIdata.lissajousPlotHandle = plot(randn(100,1));
grid off;axis off;
%-------
axes(myGUIdata.vocalTractMonitorAxis);
xdata = [0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11];
ydata = [5 5 4 4 3 3 2 2 1 1 0.2 0.2 -1 -1 -2 -2 -3 -3 -4 -4 -5 -5]/1.7;
myGUIdata.originalVTAHandle = plot(xdata,ydata*0.8+0.2,'linewidth',1,'clipping','off');
hold on
myGUIdata.manipulatedVTAHandle = plot(xdata,ydata,'linewidth',4,'color',[0 0.7 0],'clipping','off');
plot([0  11],[0 0],'clipping','off')
myGUIdata.vtlSlideHandle = plot(7.5+[-2.5 3],[2.7 2.7],'clipping','off');
myGUIdata.vtlReference = 7.5;
myGUIdata.vtlHandle = plot(myGUIdata.vtlReference,2.7,'s','linewidth',3,'markersize',10,'clipping','off');
myGUIdata.vtlMagnifier = 1;
myGUIdata.vtlFocusHandle = plot(myGUIdata.vtlReference,2.7,'o','markersize',20,'linewidth',3,'clipping','off');
set(gca,'xlim',[-0.5 10.5],'ylim',[-3 3]);
grid off;axis off;
hold off
set(myGUIdata.vtlFocusHandle,'visible','off');
set(myGUIdata.vtlText,'visible','off');
%-------
axes(myGUIdata.compositeManipulatorAxis);
xdata = (0:0.1:10)'/10;
ydata = cos(xdata*pi)+sin(xdata*pi)+sin(xdata*2*pi)+cos(xdata*pi*2);
myGUIdata.compositePlotHandle = plot(xdata,ydata,'linewidth',4,'clipping','off');
myGUIdata.rawCompositeModifier = ydata;
hold on
plot([-0.5 10.5]/10,[0 0])
knobXdata = (0:10)/10;
knobYdata = interp1(xdata,ydata,knobXdata);
knobHandle = struct;
for ii = 1:length(knobXdata)
    knobHandle(ii).handleAncor = plot(knobXdata(ii),knobYdata(ii),'sr','markersize',8,'linewidth',3);
end;
myGUIdata.knobXdata = knobXdata;
myGUIdata.knobYdata = knobYdata;
myGUIdata.knobHandle = knobHandle;
myGUIdata.knobCoefficient = zeros(length(knobXdata),1);
myGUIdata.knobFocusHandle = plot(knobXdata(1),knobYdata(1),'or','markersize',20,'linewidth',3);
myGUIdata.focusedKnobIndex = 1;
shapeBasis = zeros(length(xdata),length(knobXdata));
for ii = 1:length(knobXdata)
    shapeBasis(:,ii) = 0.5+0.5*cos(2*pi*(xdata-knobXdata(ii))*5);
    shapeBasis(abs(xdata-knobXdata(ii))>1/10,ii) = 0;
end;

%debug
if isequal(myGUIdata.debugmode,'on')
    shapeBasis 
end


myGUIdata.shapeBasis = shapeBasis;
myGUIdata.rawSectionModifier = xdata*0;
myGUIdata.shapeMagnifierSlideHandle = plot(-0.1*[1 1],[-3 3],'clipping','off');
myGUIdata.shapeMagnifierReference = 1.5;
myGUIdata.shapeMagnifierHandle = plot(-0.1,myGUIdata.shapeMagnifierReference,'s','linewidth',3,'markersize',10,'clipping','off');
myGUIdata.shapeMagnifier = 1;
myGUIdata.shapeMagnifierFocusHandle = plot(-0.1,myGUIdata.shapeMagnifierReference,'o','markersize',20,'linewidth',3,'clipping','off');
set(gca,'xlim',[-0.5 10.5]/10,'ylim',[-3 3]);
grid off;axis off;
set(myGUIdata.shapeMagnifierFocusHandle,'visible','off');
set(myGUIdata.knobFocusHandle,'visible','off');
%-------
axes(myGUIdata.individualManipulatorAxis);
xdata = (0:0.1:10)'/10;
myGUIdata.handleColor = [0.9 0 0;0 0.7 0;0 0 1;[255 62 243]/255;[51 243 144]/255;[255 155 0]/255];
myGUIdata.component1Handle = plot(xdata,cos(xdata*pi),'linewidth',3,'color',myGUIdata.handleColor(1,:),'clipping','off');
hold on
myGUIdata.component1ControlHandle = plot(0,1,'s','linewidth',3,'color',myGUIdata.handleColor(1,:),'clipping','off');
myGUIdata.component2Handle = plot(xdata,sin(xdata*pi),'linewidth',3,'color',myGUIdata.handleColor(2,:),'clipping','off');
myGUIdata.component2ControlHandle = plot(1/2,1,'s','linewidth',3,'color',myGUIdata.handleColor(2,:),'clipping','off');
myGUIdata.component3Handle = plot(xdata,sin(xdata*pi*2),'linewidth',3,'color',myGUIdata.handleColor(3,:),'clipping','off');
myGUIdata.component3ControlHandle = plot(1/4,1,'s','linewidth',3,'color',myGUIdata.handleColor(3,:),'clipping','off');
myGUIdata.component4Handle = plot(xdata,cos(xdata*pi*2),'linewidth',3,'color',myGUIdata.handleColor(4,:),'clipping','off');
myGUIdata.component4ControlHandle = plot(1,1,'s','linewidth',3,'color',myGUIdata.handleColor(4,:),'clipping','off');

myGUIdata.component5Handle = plot(xdata,sin(xdata*pi*3),'linewidth',3,'color',myGUIdata.handleColor(5,:),'clipping','off');
myGUIdata.component5ControlHandle = plot(5/6,1,'s','linewidth',3,'color',myGUIdata.handleColor(5,:),'clipping','off');
myGUIdata.component6Handle = plot(xdata,cos(xdata*pi*3),'linewidth',3,'color',myGUIdata.handleColor(6,:),'clipping','off');
myGUIdata.component6ControlHandle = plot(2/3,1,'s','linewidth',3,'color',myGUIdata.handleColor(6,:),'clipping','off');

myGUIdata.parameterFocusHandle = plot(0,1,'o','markersize',20,'linewidth',3,'clipping','off');
myGUIdata.controlHandleBundle(1).componentControlHandle = myGUIdata.component1ControlHandle;
myGUIdata.controlHandleBundle(1).componentHandle = myGUIdata.component1Handle;
myGUIdata.controlHandleBundle(2).componentControlHandle = myGUIdata.component2ControlHandle;
myGUIdata.controlHandleBundle(2).componentHandle = myGUIdata.component2Handle;

myGUIdata.controlHandleBundle(3).componentControlHandle = myGUIdata.component3ControlHandle;
myGUIdata.controlHandleBundle(3).componentHandle = myGUIdata.component3Handle;
myGUIdata.controlHandleBundle(4).componentControlHandle = myGUIdata.component4ControlHandle;
myGUIdata.controlHandleBundle(4).componentHandle = myGUIdata.component4Handle;

myGUIdata.controlHandleBundle(5).componentControlHandle = myGUIdata.component5ControlHandle;
myGUIdata.controlHandleBundle(5).componentHandle = myGUIdata.component5Handle;
myGUIdata.controlHandleBundle(6).componentControlHandle = myGUIdata.component6ControlHandle;
myGUIdata.controlHandleBundle(6).componentHandle = myGUIdata.component6Handle;

myGUIdata.numberOfBasisFunctions = 6;

plot([-0.5 10.5]/10,[0 0])
set(gca,'xlim',[-0.5 10.5]/10,'ylim',[-3 3]);
grid off;axis off;
set(myGUIdata.parameterFocusHandle,'visible','off');
%------
axes(myGUIdata.spectrogramAxis);
tx = (0:0.005:1)';
fx = (0:1024)/2048*fs;
myGUIdata.sgramImage = imagesc([tx(1) tx(end)],[0 fs/2],randn(length(fx),length(tx)));
hold on;
myGUIdata.sgramCursorHandle = plot(tx(1)*[1 1],[0 fs/2],'g','linewidth',4);
hold off;
axis('xy');
ylabel('frequency (kHz)');
hold on;
%myGUIdata.F0Axis = axes('Position',get(myGUIdata.spectrogramAxis,'Position'),...
%    'XAxisLocation','top', ...
%    'YAxisLocation','right','color','none','Xcolor','k','Ycolor','k');
%myGUIdata.F0Handle = plot(tx,tx*0+100,'color','k','parent',myGUIdata.F0Axis,'linewidth',4);
myGUIdata.F0Handle = plot(tx,tx*0+1000,'color','k','linewidth',4);
%set(myGUIdata.F0Axis,'ylim',[0 400],'xticklabel',[]);
hold off;
%get(myGUIdata.F0Axis)
set(myGUIdata.spectrogramAxis,'visible','off');
set(myGUIdata.sgramImage,'visible','off');
set(myGUIdata.sgramCursorHandle,'visible','off');
set(myGUIdata.F0Handle,'visible','off');
%set(myGUIdata.F0Axis,'visible','off');
%-------
axes(myGUIdata.tract3DAxis);
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
myGUIdata.tract3D = surf(Z,Y,X);
view(-26,12);
axis([0 1 -1 1 -1 1]);
axis off;
axis('vis3d');
%rotate3d on;
set(myGUIdata.tract3D,'visible','off');
%-------
axes(myGUIdata.spectrumMonitor);
fftl = 2048;
fx = (0:fftl-1)/fftl*4000;
myGUIdata.powerSpectrumPlot = plot(fx/1000,randn(fftl,1));
hold all
myGUIdata.originalLPCspectrumPlot = plot(fx,randn(fftl,1));
myGUIdata.modifiedLPCspectrumPlot = plot(fx,randn(fftl,1),'linewidth',2);
grid on;
xlabel('frequency (kHz)');
ylabel('level (dB)');
set(myGUIdata.spectrumMonitor,'visible','off');
set(myGUIdata.powerSpectrumPlot,'visible','off');
set(myGUIdata.originalLPCspectrumPlot,'visible','off');
set(myGUIdata.modifiedLPCspectrumPlot,'visible','off');
%-------
myGUIdata.player = audioplayer(randn(1000,1),fs);
myGUIdata.recorder = audiorecorder(fs,16,1);
myGUIdata.frameLengthInMs = 25;
myGUIdata.frameShiftInMs = 5;
%-------
set(myGUIdata.sgramDisplaySelectPanel,'visible','off');
set(myGUIdata.lissajousPlotHandle,'visible','off');
set(myGUIdata.counterText,'visible','off');
set(myGUIdata.playOriginalButton,'enable','off');
set(myGUIdata.saveOriginalButton,'enable','off');
set(myGUIdata.playSynthesizedButton,'enable','off');
set(myGUIdata.playModifiedButton,'enable','off');
set(myGUIdata.saveModifiedButton,'enable','off');
set(myGUIdata.resetMagnifierButton,'visible','off');
set(myGUIdata.f0HandleResetButton,'visible','off');
set(myGUIdata.resetVTLButton,'visible','off');
set(myGUIdata.f0RatioText,'visible','off');
set(myGUIdata.reverseValButton,'enable','off');
set(myGUIdata.loadmatfile_button,'enable','on');
myGUIdata.focus = 'off';


end

function myGUIdata = updateDisplay(myGUIdata)
set(myGUIdata.VTmainpulatorGUI,'Pointer','watch');
drawnow
%------ spectrogram
sgramStr = stftSpectrogramStructure(myGUIdata.audioData,myGUIdata.samplingFrequency,...
    myGUIdata.frameLengthInMs,myGUIdata.frameShiftInMs,'blackman');
set(myGUIdata.sgramImage,'xdata',sgramStr.temporalPositions,'ydata',sgramStr.frequencyAxis/1000, ...
    'cdata',max(-80,sgramStr.dBspectrogram));
set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
    [sgramStr.temporalPositions(1) sgramStr.temporalPositions(end)]);
set(myGUIdata.spectrogramAxis,'visible','on');
set(myGUIdata.sgramImage,'visible','on');
%------ F0
f0Struct = higherSymKalmanWithTIFupdate(myGUIdata.audioData,myGUIdata.samplingFrequency);
medianF0 = median(f0Struct.F0(f0Struct.latentSDcoeff<1.07));
myGUIdata.F0ForDisplay = f0Struct.F0/medianF0*2;
myGUIdata.f0ModRatio = 1;
set(myGUIdata.F0Handle,'visible','on','xdata',f0Struct.temporalPositions,'ydata',myGUIdata.F0ForDisplay);
set(myGUIdata.f0HandleResetButton,'visible','on');
set(myGUIdata.f0RatioText,'visible','on');
%------ information for waveform monitor (lissajous handle)
lowpassLength = round(1.7*myGUIdata.samplingFrequency/medianF0/2)*2+1;
blackmanLowPass = blackman(lowpassLength);
xxLPF = fftfilt(diff(blackmanLowPass/sum(blackmanLowPass))*100,[myGUIdata.audioData;zeros(lowpassLength,1)]);
myGUIdata.fundamentalComponent = xxLPF((1:length(myGUIdata.audioData))+floor(lowpassLength/2));
%------ LPC
convertedStr = convertTo8kHzSampledSignal(myGUIdata.audioData,myGUIdata.samplingFrequency);
myGUIdata.lpcStructure = plainLPCAnalysis(convertedStr.signal,convertedStr.samplingFrequency,myGUIdata.frameShiftInMs);
set(myGUIdata.resetVTLButton,'visible','on');
set(myGUIdata.vtlText,'visible','on');
%------ housekeeping
myGUIdata.f0Struct = f0Struct;
myGUIdata.sgramStr = sgramStr;
myGUIdata.druggingMode = 'off';
set(myGUIdata.VTmainpulatorGUI,'Pointer','arrow');
set(myGUIdata.resetMagnifierButton,'visible','on');
%------
myGUIdata = displayAreaFunctions(myGUIdata,0);
%------
set(myGUIdata.VTmainpulatorGUI,'WindowButtonMotionFcn',@WindowButtonMotionServer, ...
    'WindowButtonDownFcn',@WindowButtonDownServer,'WindowButtonUpFcn',@WindowButtonUpServer);
set(myGUIdata.sgramDisplaySelectPanel,'userdata',myGUIdata,'SelectionChangeFcn',@modeRadioButterServer,'visible','on');
end

function modeRadioButterServer(obj, event, string_arg)
myGUIdata = get(obj,'userdata');
if get(myGUIdata.stFourierTransformButton,'value') == 1
    set(myGUIdata.sgramImage,'xdata',myGUIdata.sgramStr.temporalPositions,'ydata',myGUIdata.sgramStr.frequencyAxis/1000, ...
        'cdata',max(-80,myGUIdata.sgramStr.dBspectrogram));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.sgramStr.temporalPositions(1) myGUIdata.sgramStr.temporalPositions(end)]);
elseif get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
set(obj,'userdata',myGUIdata);
end

function myGUIdata = displayModifiedLPCspectrum(myGUIdata)
responseLengthInMs = 20;
fs = 8000;
lpcStructure = myGUIdata.lpcStructure;
f0Struct = myGUIdata.f0Struct;
logArea = myGUIdata.lpcStructure.logAreaMatrix(1,:);
nSection = length(logArea);
locationList = (1:nSection)-1;
modifyerLocation = get(myGUIdata.compositePlotHandle,'xdata');
modifyerValue = get(myGUIdata.compositePlotHandle,'ydata');
sectionModifier = interp1(modifyerLocation*(nSection-1),modifyerValue,locationList,'linear','extrap');
%sgramStr = myGUIdata.sgramStr;
signalTime = (0:1/fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);
for ii = 1:fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/fs*2*pi)).*(f0InSignalTime*ii<fs/2);
end;
%outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/fs*2*pi;
fftl = length(theta);
%baseIndex = (-round(lpcStructure.frameShiftInMs/1000*fs):round(lpcStructure.frameShiftInMs/1000*fs))';
%responseLength = round(responseLengthInMs/1000*fs);
%fftLSynth = 2^ceil(log2(responseLength+length(baseIndex)+1));
%olaIndex = (1:fftLSynth)';
frequencyAxis = (0:fftl-1)'/fftl*8000;
fixedLPCspectrumInDB = zeros(fftl,length(lpcStructure.temporalPosition));
for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-0.925*cos(theta(:))) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,1)*cos(theta(:)) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,2)*cos(2*theta(:));
    logArea = lpcStructure.logAreaMatrix(ii,:);
    newLogArea = logArea(:)+sectionModifier(:)/4;
    newRef = area2ref(exp(newLogArea));
    newalp = k2alp(newRef);
    LPCspectrumLn = -2*log(abs(fft(newalp(:),fftl)));
    %LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedLPCspectrumInPower = exp(fixedLPCspectrum);
    fixedLPCspectrumInPower = fixedLPCspectrumInPower ...
        /sum(fixedLPCspectrumInPower)*lpcStructure.powerList(ii);% /f0Struct.F0(ii); for display, this is not needed.
    fixedLPCspectrumInDB(:,ii) = 10*log10(fixedLPCspectrumInPower(:));
end;
fixedLPCspectrumInDB = fixedLPCspectrumInDB-max(fixedLPCspectrumInDB(:));
fixedLPCspectrumInDB = interp1(frequencyAxis,fixedLPCspectrumInDB, ...
    frequencyAxis*myGUIdata.vtlMagnifier,'linear',-100);
myGUIdata.fixedLPCspectrumInDB = fixedLPCspectrumInDB(1:fftl/2+1,:)-max(fixedLPCspectrumInDB(:));
end

function straightSgram = sgram2STRAIGHTSgram(f0Struct,sgramStr)
%tx = sgramStr
straightSgram.sgram = sgramStr;
straightSgram.f0 = f0Struct;
end

% ---
function WindowButtonMotionServer(src, evnt)
myGUIdata = guidata(src);
currentPointInParam = get(myGUIdata.individualManipulatorAxis,'CurrentPoint');
currentPointInSgram = get(myGUIdata.spectrogramAxis,'CurrentPoint');
currentPointInShape = get(myGUIdata.compositeManipulatorAxis,'CurrentPoint');
currentPointInLogArea = get(myGUIdata.vocalTractMonitorAxis,'CurrentPoint');
xlim = get(myGUIdata.individualManipulatorAxis,'xlim');
ylim = get(myGUIdata.individualManipulatorAxis,'ylim');
controlHandlePosition = zeros(3,2);
xlimSgram = get(myGUIdata.spectrogramAxis,'xlim');
ylimSgram = get(myGUIdata.spectrogramAxis,'ylim');
xlimShape = get(myGUIdata.compositeManipulatorAxis,'xlim');
xlimShape(1) = -0.15;
ylimShape = get(myGUIdata.compositeManipulatorAxis,'ylim');
xlimLogArea = get(myGUIdata.vocalTractMonitorAxis,'xlim');
ylimLogArea = get(myGUIdata.vocalTractMonitorAxis,'ylim');
positionHolder = zeros(2,myGUIdata.numberOfBasisFunctions);
switch myGUIdata.druggingMode
    case 'off'
        xinit = currentPointInParam(1,1);
        yinit = currentPointInParam(1,2);
        positionHolder(1,:) = xinit;
        positionHolder(2,:) = yinit;
        xinitSgram = currentPointInSgram(1,1);
        yinitSgram = currentPointInSgram(1,2);
        xinitShape = currentPointInShape(1,1);
        yinitShape = currentPointInShape(1,2);
        xinitLogArea = currentPointInLogArea(1,1);
        yinitLogArea = currentPointInLogArea(1,2);
        
        %indivisual plain
        if (xlim(1)-xinit)*(xlim(2)-xinit) <= 0 && (ylim(1)-yinit)*(ylim(2)-yinit) <= 0
            set(myGUIdata.VTmainpulatorGUI,'Pointer','fullcrosshair');
            %set(myGUIdata.parameterFocusHandle,'xdata',xinit,'ydata',yinit);
            for ii = 1:myGUIdata.numberOfBasisFunctions
                controlHandlePosition(ii,1) = get(myGUIdata.controlHandleBundle(ii).componentControlHandle,'xdata');
                controlHandlePosition(ii,2) = get(myGUIdata.controlHandleBundle(ii).componentControlHandle,'ydata');
            end;
            distanceVector = sqrt(sum((controlHandlePosition'-positionHolder).^2));
            [minimumDist,minimumIndex] = min(distanceVector);
            if minimumDist < 0.08
                xFocus = controlHandlePosition(minimumIndex,1);
                yFocus = controlHandlePosition(minimumIndex,2);
                set(myGUIdata.parameterFocusHandle,'xdata',xFocus,'ydata',yFocus, ...
                    'visible','on','color',myGUIdata.handleColor(minimumIndex,:));
                for ii = 1:myGUIdata.numberOfBasisFunctions
                    if ii == minimumIndex
                        currentWidth = 6;
                    else
                        currentWidth = 2;
                    end
                    set(myGUIdata.controlHandleBundle(ii).componentHandle,'linewidth',currentWidth);
                end;
                myGUIdata.focus = 'on';
                myGUIdata.xFocus = xFocus;
                myGUIdata.focusIndex = minimumIndex;
            else
                set(myGUIdata.parameterFocusHandle,'visible','off');
                for ii = 1:myGUIdata.numberOfBasisFunctions
                    set(myGUIdata.controlHandleBundle(ii).componentHandle,'linewidth',2);
                end;
            end;
        %spectrogram plain
        elseif (xlimSgram(1)-xinitSgram)*(xlimSgram(2)-xinitSgram) <= 0 && ...
                (ylimSgram(1)-yinitSgram)*(ylimSgram(2)-yinitSgram) <= 0
            xF0data = get(myGUIdata.F0Handle,'xdata');
            yF0data = get(myGUIdata.F0Handle,'ydata');
            [minXdata,indexXdata] = min(abs(xF0data-xinitSgram));
            if abs(yF0data(indexXdata)-yinitSgram) < 0.15
                set(myGUIdata.VTmainpulatorGUI,'Pointer','hand');
            else
                set(myGUIdata.VTmainpulatorGUI,'Pointer','cross');
                myGUIdata.focus = 'on';
            end
        %composite plain
        elseif  (xlimShape(1)-xinitShape)*(xlimShape(2)-xinitShape) <= 0 && ...
                (ylimShape(1)-yinitShape)*(ylimShape(2)-yinitShape) <= 0
            yShapeKnobData = get(myGUIdata.shapeMagnifierHandle,'ydata');
            %set(myGUIdata.VTmainpulatorGUI,'Pointer','circle');drawnow
            if abs(yShapeKnobData-yinitShape) < 0.15 && ...
                    abs(xinitShape-(-0.1)) < 0.015
                set(myGUIdata.VTmainpulatorGUI,'Pointer','crosshair');
                set(myGUIdata.shapeMagnifierFocusHandle,'visible','on');
                myGUIdata.focus = 'on';
                %get(myGUIdata.VTmainpulatorGUI,'Pointer')
            elseif min(abs(myGUIdata.knobXdata-xinitShape)) < 0.015
                [~,knobIndex] = min(abs(myGUIdata.knobXdata-xinitShape));
                if abs(myGUIdata.knobYdata(knobIndex)-yinitShape) < 0.15
                    customCursor(myGUIdata.VTmainpulatorGUI);
                    set(myGUIdata.knobFocusHandle,'visible','on','xdata',myGUIdata.knobXdata(knobIndex), ...
                        'ydata',myGUIdata.knobYdata(knobIndex));
                    myGUIdata.focus = 'on';
                    myGUIdata.focusedKnobIndex = knobIndex;
                else
                    myGUIdata.focus = 'off';
                    set(myGUIdata.VTmainpulatorGUI,'Pointer','arrow');
                    set(myGUIdata.knobFocusHandle,'visible','off');
                end;
            else
                set(myGUIdata.shapeMagnifierFocusHandle,'visible','off');
                set(myGUIdata.VTmainpulatorGUI,'Pointer','arrow');
                myGUIdata.focus = 'off';
                %drawnow
            end;
        %vtl plain
        elseif  (xlimLogArea(1)-xinitLogArea)*(xlimLogArea(2)-xinitLogArea) <= 0 && ...
                (ylimLogArea(1)-yinitLogArea)*(ylimLogArea(2)-yinitLogArea) <= 0
            xVTLKnobData = get(myGUIdata.vtlHandle,'xdata');
            %set(myGUIdata.VTmainpulatorGUI,'Pointer','circle');drawnow
            if abs(xVTLKnobData-xinitLogArea) < 0.15 && ...
                    abs(yinitLogArea-(2.7)) < 0.15
                set(myGUIdata.VTmainpulatorGUI,'Pointer','circle');
                set(myGUIdata.vtlFocusHandle,'visible','on');
                myGUIdata.focus = 'on';
                %get(myGUIdata.VTmainpulatorGUI,'Pointer')
            else
                set(myGUIdata.vtlFocusHandle,'visible','off');
                set(myGUIdata.VTmainpulatorGUI,'Pointer','arrow');
                myGUIdata.focus = 'off';
                %drawnow
            end;
        else
            set(myGUIdata.VTmainpulatorGUI,'Pointer','arrow');
        end
        guidata(src,myGUIdata);
    case 'on'
        %xinit = currentPointInParam(1,1);
        switch get(myGUIdata.VTmainpulatorGUI,'Pointer');
            case 'fullcrosshair'    %indivisual
                yinit = currentPointInParam(1,2);
                xdata = get(myGUIdata.controlHandleBundle(myGUIdata.focusIndex).componentHandle,'xdata');
                set(myGUIdata.controlHandleBundle(myGUIdata.focusIndex).componentControlHandle,'ydata',yinit);
                set(myGUIdata.controlHandleBundle(myGUIdata.focusIndex).componentControlHandle,'xdata',myGUIdata.xFocus);
                switch myGUIdata.focusIndex
                    case 1
                        ydata = yinit*cos(xdata*pi);
                    case 2
                        ydata = yinit*sin(xdata*pi);
                    case 3
                        ydata = yinit*sin(2*xdata*pi);
                    case 4
                        ydata = yinit*cos(2*xdata*pi);
                    case 5
                        ydata = yinit*sin(3*xdata*pi);
                    case 6
                        ydata = yinit*cos(3*xdata*pi);
                end
                set(myGUIdata.controlHandleBundle(myGUIdata.focusIndex).componentHandle,'ydata',ydata);
                set(myGUIdata.parameterFocusHandle,'ydata',yinit);
                for ii = 1:myGUIdata.numberOfBasisFunctions
                    if ii ~= myGUIdata.focusIndex
                        ydata = ydata+get(myGUIdata.controlHandleBundle(ii).componentHandle,'ydata');
                    end;
                end;
                
                myGUIdata.rawCompositeModifier = ydata(:);
                set(myGUIdata.compositePlotHandle,'ydata',...
                    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
                myGUIdata = updateKnobs(myGUIdata);
                myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
                
                if get(myGUIdata.autoRegressiveModelButton,'value') == 1
                    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
                    fftl = length(myGUIdata.lpcStructure.rawFAxis);
                    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
                        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
                        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
                    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
                        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
                end;
                
                guidata(src,myGUIdata);
            case 'cross'    %spectrogram
                xinit = currentPointInSgram(1,1);
                set(myGUIdata.sgramCursorHandle,'visible','on','xdata', ...
                    max(xlimSgram(1),min(xlimSgram(2),xinit))*[1 1]);
                myGUIdata = displayAreaFunctions(myGUIdata,xinit);
                guidata(src,myGUIdata);
            case 'hand' %F0
                xF0data = get(myGUIdata.F0Handle,'xdata');
                [minXdata,indexXdata] = min(abs(xF0data-currentPointInSgram(1,1)));
                myGUIdata.f0ModRatio = max(0.167,min(6,currentPointInSgram(1,2)/myGUIdata.F0ForDisplay(indexXdata)));
                set(myGUIdata.F0Handle,'ydata',myGUIdata.f0ModRatio*myGUIdata.F0ForDisplay);
                set(myGUIdata.f0RatioText,'string',num2str(myGUIdata.f0ModRatio,'%4.2f'),'visible','on');
                guidata(src,myGUIdata);
            case 'custom'   %composite Knob
                %myGUIdata.focusedKnobIndex
                yinitShape = currentPointInShape(1,2);
                set(myGUIdata.knobFocusHandle,'visible','on','ydata',yinitShape);
                set(myGUIdata.knobHandle(myGUIdata.focusedKnobIndex).handleAncor,'ydata',yinitShape);
                xdata = get(myGUIdata.compositePlotHandle,'xdata');
                myGUIdata.knobYdata(myGUIdata.focusedKnobIndex) = yinitShape;
                
                %debug
                if isequal(myGUIdata.debugmode,'on')
                    disp('rawCompositeModifier')
                    myGUIdata.rawCompositeModifier(1:10)
                    disp('rawSection')
                    myGUIdata.rawSectionModifier(1:10)
                end
                
                tmpVisibleModifierValue = interp1(xdata,myGUIdata.rawCompositeModifier*myGUIdata.shapeMagnifier, ...
                    myGUIdata.knobXdata(myGUIdata.focusedKnobIndex));
                if myGUIdata.shapeMagnifier ~= 0
                    myGUIdata.knobCoefficient(myGUIdata.focusedKnobIndex) = ...
                        (yinitShape-tmpVisibleModifierValue)/myGUIdata.shapeMagnifier;
                end;
                
                myGUIdata = calculateSectionModifier(myGUIdata);
                set(myGUIdata.compositePlotHandle,'ydata',...
                    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
                
                myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
                
                if get(myGUIdata.autoRegressiveModelButton,'value') == 1
                    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
                    fftl = length(myGUIdata.lpcStructure.rawFAxis);
                    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
                        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
                        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
                    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
                        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
                end;
                guidata(src,myGUIdata);
            case 'crosshair'    %composite 
                myGUIdata.shapeMagnifier = max(ylimShape(1),min(ylimShape(2), ...
                    currentPointInShape(1,2)))/myGUIdata.shapeMagnifierReference;
                set(myGUIdata.shapeMagnifierHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
                set(myGUIdata.shapeMagnifierFocusHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
                set(myGUIdata.shapeMagnifierText,'string',num2str(myGUIdata.shapeMagnifier,'%4.2f'));
                set(myGUIdata.compositePlotHandle,'ydata',...
                    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
                myGUIdata = updateKnobs(myGUIdata);
                myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
                if get(myGUIdata.autoRegressiveModelButton,'value') == 1
                    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
                    fftl = length(myGUIdata.lpcStructure.rawFAxis);
                    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
                        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
                        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
                    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
                        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
                end;
                guidata(src,myGUIdata);
            case 'circle'   %vtl
                myGUIdata.vtlMagnifier = max(5,min(10.5, ...
                    currentPointInLogArea(1,1)))/myGUIdata.vtlReference;
                set(myGUIdata.vtlHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
                set(myGUIdata.vtlFocusHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
                set(myGUIdata.vtlText,'string',num2str(myGUIdata.vtlMagnifier,'%4.2f'));
                xdata = get(myGUIdata.originalVTAHandle,'xdata');
                set(myGUIdata.manipulatedVTAHandle,'xdata',xdata*myGUIdata.vtlMagnifier);
                myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
                if 1 == 1
                    if get(myGUIdata.autoRegressiveModelButton,'value') == 1
                        myGUIdata = displayModifiedLPCspectrum(myGUIdata);
                        fftl = length(myGUIdata.lpcStructure.rawFAxis);
                        set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
                            myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
                            'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
                        set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
                            [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
                    end;
                end;
                guidata(src,myGUIdata);
        end;
end;
end

function myGUIdata = updateKnobs(myGUIdata)
xdata = get(myGUIdata.compositePlotHandle,'xdata');
ydata = get(myGUIdata.compositePlotHandle,'ydata');
nKnobs = length(myGUIdata.knobCoefficient);
for ii = 1:nKnobs
    myGUIdata.knobYdata(ii) = interp1(xdata,ydata,myGUIdata.knobXdata(ii));
    set(myGUIdata.knobHandle(ii).handleAncor,'ydata',myGUIdata.knobYdata(ii));
end
end

function myGUIdata = calculateSectionModifier(myGUIdata)
nSection = length(myGUIdata.knobCoefficient);
rawSectionModifier = myGUIdata.rawCompositeModifier(:)*0;

%debug
if isequal(myGUIdata.debugmode,'on')
    myGUIdata.shapeBasis
    
end

for ii = 1:nSection
    rawSectionModifier = rawSectionModifier+myGUIdata.knobCoefficient(ii)*myGUIdata.shapeBasis(:,ii);
end;
myGUIdata.rawSectionModifier = rawSectionModifier;
end

function myGUIdata = displayAreaFunctions(myGUIdata,xinit)
[minValue,minIndex] = min(abs(myGUIdata.lpcStructure.temporalPosition-xinit));
[minSgramValue,minSgramIndex] = min(abs(myGUIdata.sgramStr.temporalPositions-xinit));
logArea = myGUIdata.lpcStructure.logAreaMatrix(minIndex,:);
nSection = length(logArea);
locationList = (1:nSection)-1;
xdata = get(myGUIdata.originalVTAHandle,'xdata');
displayLogArea = interp1(locationList,logArea,xdata,'nearest','extrap');
displayLogArea = displayLogArea(:);
displayLogArea = [displayLogArea(1); displayLogArea(1:end-1)];
set(myGUIdata.originalVTAHandle,'ydata',displayLogArea-mean(displayLogArea));
modifyerLocation = get(myGUIdata.compositePlotHandle,'xdata');
modifyerValue = get(myGUIdata.compositePlotHandle,'ydata');
sectionModifier = interp1(modifyerLocation*(nSection-2),modifyerValue,locationList,'linear','extrap');
myGUIdata.sectionModifier = sectionModifier;
displayLogArea = interp1(locationList,logArea(:)+sectionModifier(:)/4,xdata,'nearest','extrap');
displayLogArea = displayLogArea(:);
displayLogArea = [displayLogArea(1); displayLogArea(1:end-1)];
set(myGUIdata.manipulatedVTAHandle,'ydata',displayLogArea-mean(displayLogArea));
[X,Y,Z] = cylinder(tubeDisplay(logArea(:)+sectionModifier(:)/4),40);
set(myGUIdata.tract3D,'xdata',Z,'ydata',Y,'zdata',X,'visible','on');
newLogArea = logArea(:)+sectionModifier(:)/4;
newRef = area2ref(exp(newLogArea));
newalp = k2alp(newRef);
newLPCspectrumLn = -2*log(abs(fft(newalp,myGUIdata.lpcStructure.fftl)));
theta = myGUIdata.lpcStructure.rawFAxis/8000*2*pi;
%inverseHPF = -2*log(1-0.925*cos(theta(:)));
inversePreprocessingShape = -2*log(1-myGUIdata.lpcStructure.preEmphasis*cos(theta(:))) ...
    +myGUIdata.lpcStructure.cepstrumList(minIndex,1)*cos(theta(:)) ...
    +myGUIdata.lpcStructure.cepstrumList(minIndex,2)*cos(2*theta(:));
%fixedLPCspectrum = newLPCspectrumLn(:)+myGUIdata.lpcStructure.cepstrumList(minIndex,1)*cos(theta(:)) ...
%    +myGUIdata.lpcStructure.cepstrumList(minIndex,2)*cos(2*theta(:))+inverseHPF;
%fixedLPCspectrumInDB = 10*fixedLPCspectrum/log(10);
fixedLPCspectrum = newLPCspectrumLn(:)+inversePreprocessingShape;
fixedLPCspectrumInPower = exp(fixedLPCspectrum);
dBspectrumInPower = 10.0.^(myGUIdata.sgramStr.dBspectrogram(:,minSgramIndex)/10);
fixedLPCspectrumInPower = fixedLPCspectrumInPower/mean(fixedLPCspectrumInPower) ...
    *mean(dBspectrumInPower(myGUIdata.sgramStr.frequencyAxis<4000));
fixedLPCspectrumInDB = 10*log10(fixedLPCspectrumInPower);
%fixedLPCspectrumInDB = fixedLPCspectrumInDB-mean(fixedLPCspectrumInDB)+mean(myGUIdata.sgramStr.dBspectrogram(:,minSgramIndex));
set(myGUIdata.spectrumMonitor,'visible','on','ylim',[-80 10],'xlim',[0 4000]);
set(myGUIdata.powerSpectrumPlot,'visible','on','xdata', ...
    myGUIdata.sgramStr.frequencyAxis,'ydata',myGUIdata.sgramStr.dBspectrogram(:,minSgramIndex));
set(myGUIdata.originalLPCspectrumPlot,'visible','off');
set(myGUIdata.modifiedLPCspectrumPlot,'visible','on','xdata',myGUIdata.lpcStructure.rawFAxis/myGUIdata.vtlMagnifier,...
    'ydata',fixedLPCspectrumInDB);
myGUIdata.xinit = xinit;
end

function crossSection = tubeDisplay(logArea)
areaList = exp(logArea);
crossSection = sqrt(areaList/sum(areaList));
end

function WindowButtonDownServer(src, evnt)
myGUIdata = guidata(src);
switch myGUIdata.focus
    case 'on'
        myGUIdata.druggingMode = 'on';
        guidata(src,myGUIdata);
end;
switch get(myGUIdata.VTmainpulatorGUI,'Pointer')
    case 'hand'
        set(myGUIdata.F0Handle,'color',[0.3 1 0])
    case 'crosshair'
        set(myGUIdata.shapeMagnifierSlideHandle,'linewidth',5,'color',[0.3 1 0]);
    case 'circle'
        set(myGUIdata.vtlSlideHandle,'linewidth',5,'color',[0.3 1 0]);
    case 'custom'
        %set(myGUIdata.knobFocusHandle,'color',[0.3 1 0]);
        set(myGUIdata.knobFocusHandle,'color',[1 0.1 0.1]);
        %set(myGUIdata.knobHandle(myGUIdata.focusedKnobIndex).handleAncor,'color',[0.3 1 0]);
        set(myGUIdata.knobHandle(myGUIdata.focusedKnobIndex).handleAncor,'color',[0.1 0.1 1]);
end
end

function WindowButtonUpServer(src, evnt)
myGUIdata = guidata(src);
myGUIdata.druggingMode = 'off';
myGUIdata.focus = 'off';
set(myGUIdata.F0Handle,'color',[0 0 0])
set(myGUIdata.shapeMagnifierSlideHandle,'linewidth',1,'color',[0 0 1]);
set(myGUIdata.knobHandle(myGUIdata.focusedKnobIndex).handleAncor,'color',[0 0 1]);
set(myGUIdata.knobFocusHandle,'color',[0 0 1]);
set(myGUIdata.vtlSlideHandle,'linewidth',1,'color',[0 0 1]);
%set(myGUIdata.sgramCursorHandle,'visible','off');
%set(myGUIdata.tract3D,'visible','off');
guidata(src,myGUIdata);
end

% --- Executes on button press in recordButton.
function recordButton_Callback(hObject, eventdata, handles)
% hObject    handle to recordButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
set(myGUIdata.counterText,'visible','off');
set(myGUIdata.playOriginalButton,'enable','off');
set(myGUIdata.saveOriginalButton,'enable','off');
set(myGUIdata.playSynthesizedButton,'enable','off');
set(myGUIdata.playModifiedButton,'enable','off');
set(myGUIdata.saveModifiedButton,'enable','off');
myGUIdata.druggingMode = 'off';
switch get(myGUIdata.recorder,'running')
    case 'off'
        myGUIdata.recorder = audiorecorder(myGUIdata.samplingFrequency,16,1);
        myGUIdata.audoRecordTimerPeriod = 0.053;
        set(myGUIdata.recordButton,'string','STOP');
        set(myGUIdata.recorder,'userdata',myGUIdata.VTmainpulatorGUI,'TimerFcn', ...
            @monitorUpdateTimerFunction,'TimerPeriod',myGUIdata.audoRecordTimerPeriod, ...
            'StopFcn',@recordingStopFcn);
        set(myGUIdata.lissajousPlotHandle,'visible','on');
        set(myGUIdata.counterText,'visible','on');
        set(myGUIdata.counterText,'string','00000 ms');
        guidata(handles.VTmainpulatorGUI,myGUIdata);
        record(myGUIdata.recorder,5);
        guidata(handles.VTmainpulatorGUI,myGUIdata);
    case 'on'
        %myGUIdata.audioData = getaudiodata(myGUIdata.recorder);
        stop(myGUIdata.recorder);
end
end

function recordingStopFcn(obj, event, string_arg)
handleForTimer = get(obj,'userData');
myGUIdata = guidata(handleForTimer);
myGUIdata.audioData = getaudiodata(myGUIdata.recorder);
set(myGUIdata.recordButton,'string','RECORD');
set(myGUIdata.lissajousPlotHandle,'visible','off');
set(myGUIdata.counterText,'visible','off');
%----
set(myGUIdata.stFourierTransformButton,'value',1);
myGUIdata = updateDisplay(myGUIdata);
%----
myGUIdata.dataDescription = ['new recording:' datestr(now)];
set(myGUIdata.playOriginalButton,'enable','on');
set(myGUIdata.saveOriginalButton,'enable','on');
set(myGUIdata.playSynthesizedButton,'enable','on');
set(myGUIdata.playModifiedButton,'enable','on'); %TNX! Yoshimoto-san,
guidata(handleForTimer,myGUIdata);
end

function monitorUpdateTimerFunction(obj, event, string_arg)
%   timer function for updating display

%   14/May/2014
%   by Hideki Kawahara
handleForTimer = get(obj,'userData');
%myGUIdata = get(handleForTimer,'userdata');
myGUIdata = guidata(handleForTimer);
currentPoint = get(myGUIdata.recorder,'CurrentSample');
fs = get(myGUIdata.recorder,'SampleRate');
index30ms = -round(0.03*fs):0;
%TotalSamples = get(myGUIdata.recorder,'TotalSamples');
y = getaudiodata(myGUIdata.recorder);
ydata = y(max(1,min(length(y),length(y)+index30ms)));
set(myGUIdata.lissajousAxis,'ylim',max(abs(ydata))*[-1 1],'xlim',[index30ms(1) index30ms(end)]);
set(myGUIdata.lissajousPlotHandle,'xdata',index30ms,'ydata',ydata);
set(myGUIdata.counterText,'string',[num2str(currentPoint/fs*1000,'%05.0f') '  ms']);
end

% --- Executes on button press in loadFileButton.
function loadFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
[file,path] = uigetfile('*.wav','Select sound file');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        disp('Load is cancelled!');
        return;
    end;
end;

[x,fs] = audioread([path file]);

myGUIdata.audioData = x(:,1);
myGUIdata.samplingFrequency = fs;
myGUIdata.dataDescription = ['file name:' file ' path:' path];

%STRIAHGT analysis
%r1 = extractInitialF0bySymmetryICSP2013(x,fsForAnalysis);
r1 = extractInitialF0bySymmetryICSP2013(x,fs);
f0 = r1.f0Raw*0+median(r1.f0Raw(r1.rawPeriodicity>0/7))*0.7;
periodicityLevel = f0*0+0.9;
%xClean = removeLF(x,fsForAnalysis,f0,periodicityLevel);
xClean = removeLF(x,fs,f0,periodicityLevel);
%r1 = extractInitialF0bySymmetryICSP2013(xClean,fsForAnalysis);
r1 = extractInitialF0bySymmetryICSP2013(xClean,fs);
rClean = cleaningF0RawTrajectoryV2(r1);
%refR1 = refineExGxOutputRevSqx(xClean,fsForAnalysis,rClean);
refR1 = refineExGxOutputRevSqx(xClean,fs,rClean);
rClean.f0 = refR1.f0Refined;
rClean.f0(refR1.rawPeriodicity<0.3) = median(rClean.f0(refR1.rawPeriodicity>0.7));
rClean.samplingFrequency = fs;
%f1 = exSpectrumTSTRAIGHTGB(xSegment,fsForAnalysis,rClean);
myGUIdata.f1 = exSpectrumTSTRAIGHTGB(x,fs,rClean);

%{
%rClean.vuv = refineVoicingDecision(x,rClean);
q = aperiodicityRatioSigmoid(x,rClean,1,2,0); % aperiodicity extraction
f = exSpectrumTSTRAIGHTGB(x,fs,q);
STRAIGHTobject.waveform = x;
STRAIGHTobject.samplingFrequency = fs;
STRAIGHTobject.refinedF0Structure.temporalPositions = r.temporalPositions;
STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT = f.spectrogramSTRAIGHT;
STRAIGHTobject.refinedF0Structure.vuv = rClean.vuv;
f.spectrogramSTRAIGHT = unvoicedProcessing(STRAIGHTobject);
s2 = exGeneralSTRAIGHTsynthesisR2(q,f); % new implementation
sound(s2.synthesisOut/max(abs(s2.synthesisOut))*0.8,fs)
%}
%----
set(myGUIdata.stFourierTransformButton,'value',1);
myGUIdata = updateDisplay(myGUIdata);
%----
set(myGUIdata.playOriginalButton,'enable','on');
set(myGUIdata.playSynthesizedButton,'enable','on');
set(myGUIdata.saveModifiedButton,'enable','off');
set(myGUIdata.saveOriginalButton,'enable','off');
set(myGUIdata.playModifiedButton,'enable','on');
guidata(handles.VTmainpulatorGUI,myGUIdata);
end

% --- Executes on button press in loadmatfile_button.
function loadmatfile_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadmatfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
[file,path] = uigetfile('*.mat','Select sound file');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        disp('Load is cancelled!');
        return;
    end;
end;
load([path file])

myGUIdata.MatDataDescription = ['file name:' file ' path:' path];

myGUIdata.audioData = modificationStructure.originalSignal;
myGUIdata.samplingFrequency = modificationStructure.samplingFrequency;
myGUIdata.dataDescription = modificationStructure.dataDescription;
myGUIdata.lpcStructure = modificationStructure.lpcStructure;
myGUIdata.f0Struct = modificationStructure.f0Struct;
myGUIdata.sectionModifier = modificationStructure.sectionModifier;
myGUIdata.synthesizedSignal = modificationStructure.synthesizedSignal;
myGUIdata.samplingFrequencyOfSynth = modificationStructure.samplingFrequencyOfSynth;
myGUIdata.f0ModRatio = modificationStructure.f0ModificationRatio;
myGUIdata.rawCompositeModifier = modificationStructure.rawCompositeModifier;
myGUIdata.rawSectionModifier = modificationStructure.rawSectionModifier;
myGUIdata.knobCoefficient = modificationStructure.knobCoefficient;
myGUIdata.shapeMagnifier = modificationStructure.shapeMagnifier;
myGUIdata.vtlMagnifier = modificationStructure.vtlMagnifier;

myGUIdata = updateDisplay(myGUIdata);

%-----------Indivisual
%{
for focusIndex = 1:myGUIdata.numberOfBasisFunctions
    xdata = get(myGUIdata.controlHandleBundle(focusIndex).componentHandle,'xdata');
    ydata = xdata*0;
    set(myGUIdata.controlHandleBundle(focusIndex).componentControlHandle,'ydata',0);
    set(myGUIdata.controlHandleBundle(focusIndex).componentHandle,'ydata',ydata);
end;
myGUIdata.rawCompositeModifier = ydata(:);
set(myGUIdata.compositePlotHandle,'ydata',(myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
%set(myGUIdata.compositePlotHandle,'ydata',ydata);
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
myGUIdata = displayModifiedLPCspectrum(myGUIdata);
fftl = length(myGUIdata.lpcStructure.rawFAxis);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
%}
%------------composite
set(myGUIdata.shapeMagnifierHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
set(myGUIdata.shapeMagnifierFocusHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
set(myGUIdata.shapeMagnifierText,'string',num2str(myGUIdata.shapeMagnifier,'%4.2f'));
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);
set(myGUIdata.compositePlotHandle,'ydata',...
    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);

if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;

%------------VTL
set(myGUIdata.vtlHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
myGUIdata.vtlMagnifier = modificationStructure.vtlMagnifier;
set(myGUIdata.vtlText,'string',num2str(myGUIdata.vtlMagnifier,'%4.2f'));
set(myGUIdata.vtlFocusHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
xdata = get(myGUIdata.originalVTAHandle,'xdata');
set(myGUIdata.manipulatedVTAHandle,'xdata',xdata*myGUIdata.vtlMagnifier);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;


%------------F0
set(myGUIdata.F0Handle,'ydata',myGUIdata.f0ModRatio*myGUIdata.F0ForDisplay);
myGUIdata.f0ModRatio = modificationStructure.f0ModificationRatio;
set(myGUIdata.f0RatioText,'string',num2str(myGUIdata.f0ModRatio,'%4.2f'),'visible','on');

%------------Knob
myGUIdata.knobCoefficient = modificationStructure.knobCoefficient;
myGUIdata = calculateSectionModifier(myGUIdata);
set(myGUIdata.compositePlotHandle,'ydata',...
    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);

if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;

%--------------------------------
set(myGUIdata.playOriginalButton,'enable','on');
set(myGUIdata.playSynthesizedButton,'enable','on');
set(myGUIdata.saveModifiedButton,'enable','off');
set(myGUIdata.saveOriginalButton,'enable','off');
set(myGUIdata.playModifiedButton,'enable','on');
set(myGUIdata.reverseValButton,'enable','on');
guidata(handles.VTmainpulatorGUI,myGUIdata);

end


% --- Executes on button press in reverseValButton.
function reverseValButton_Callback(hObject, eventdata, handles)
% hObject    handle to reverseValButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);

%{
zeroind = find(myGUIdata.sectionModifier == 0);
nzeroind = find(myGUIdata.sectionModifier);
if ~isempty(zeroind)
    myGUIdata.sectionModifier(nzeroind) = -myGUIdata.sectionModifier(nzeroind);
else
    myGUIdata.sectionModifier = -myGUIdata.sectionModifier;
end

zeroind = find(myGUIdata.rawSectionModifier == 0);
nzeroind = find(myGUIdata.rawSectionModifier);
if ~isempty(zeroind)
    myGUIdata.rawSectionModifier(nzeroind) = -myGUIdata.rawSectionModifier(nzeroind);
else
    myGUIdata.rawSectionModifier = -myGUIdata.rawSectionModifier;
end

zeroind = find(myGUIdata.rawCompositeModifier == 0);
nzeroind = find(myGUIdata.rawCompositeModifier);
if ~isempty(zeroind)
    myGUIdata.rawCompositeModifier(nzeroind) = -myGUIdata.rawCompositeModifier(nzeroind);
else
    myGUIdata.rawCompositeModifier = -myGUIdata.rawCompositeModifier;
end

beforeknobCoefficient = myGUIdata.knobCoefficient;
zeroind = find(myGUIdata.knobCoefficient == 0);
nzeroind = find(myGUIdata.knobCoefficient);
if ~isempty(zeroind)
    %myGUIdata.knobCoefficient(nzeroind) = 1./myGUIdata.knobCoefficient(nzeroind);
    myGUIdata.knobCoefficient(nzeroind) = -myGUIdata.knobCoefficient(nzeroind);
else
    %myGUIdata.knobCoefficient = 1./myGUIdata.knobCoefficient;
    myGUIdata.knobCoefficient = -myGUIdata.knobCoefficient;
end

%}

%set matFile field

myGUIdata.sectionModifier = myGUIdata.sectionModifier;
myGUIdata.rawCompositeModifier = myGUIdata.rawCompositeModifier;
myGUIdata.knobCoefficient = myGUIdata.knobCoefficient;

if myGUIdata.f0ModRatio ~= 0 
    myGUIdata.f0ModRatio = 1/myGUIdata.f0ModRatio;
end
if myGUIdata.shapeMagnifier ~= 0
    myGUIdata.shapeMagnifier = -myGUIdata.shapeMagnifier;
end
if myGUIdata.vtlMagnifier ~= 0
    myGUIdata.vtlMagnifier = 1/myGUIdata.vtlMagnifier;
end


%------------Knob
%myGUIdata.knobCoefficient = modificationStructure.knobCoefficient;
myGUIdata = calculateSectionModifier(myGUIdata);
set(myGUIdata.compositePlotHandle,'ydata',...
    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);

myGUIdata = updateKnobs(myGUIdata);

%------------composite
set(myGUIdata.shapeMagnifierHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
set(myGUIdata.shapeMagnifierFocusHandle,'ydata',myGUIdata.shapeMagnifier*myGUIdata.shapeMagnifierReference);
set(myGUIdata.shapeMagnifierText,'string',num2str(myGUIdata.shapeMagnifier,'%4.2f'));

myGUIdata = calculateSectionModifier(myGUIdata);
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);
set(myGUIdata.compositePlotHandle,'ydata',...
    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);

%myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);


%------------VTL
set(myGUIdata.vtlHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
%myGUIdata.vtlMagnifier = modificationStructure.vtlMagnifier;
set(myGUIdata.vtlText,'string',num2str(myGUIdata.vtlMagnifier,'%4.2f'));
set(myGUIdata.vtlFocusHandle,'xdata',myGUIdata.vtlMagnifier*myGUIdata.vtlReference);
xdata = get(myGUIdata.originalVTAHandle,'xdata');
set(myGUIdata.manipulatedVTAHandle,'xdata',xdata*myGUIdata.vtlMagnifier);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);

%------------F0
set(myGUIdata.F0Handle,'ydata',myGUIdata.f0ModRatio*myGUIdata.F0ForDisplay);
%myGUIdata.f0ModRatio = modificationStructure.f0ModificationRatio;
set(myGUIdata.f0RatioText,'string',num2str(myGUIdata.f0ModRatio,'%4.2f'),'visible','on');


myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;

%--------------------------------
set(myGUIdata.playOriginalButton,'enable','on');
set(myGUIdata.playSynthesizedButton,'enable','on');
set(myGUIdata.saveModifiedButton,'enable','off');
set(myGUIdata.saveOriginalButton,'enable','off');
set(myGUIdata.playModifiedButton,'enable','on');
set(myGUIdata.reverseValButton,'enable','on');
set(myGUIdata.reverseValButton,'enable','off');
guidata(handles.VTmainpulatorGUI,myGUIdata);

end


% --- Executes on button press in playOriginalButton.
function playOriginalButton_Callback(hObject, eventdata, handles)
% hObject    handle to playOriginalButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
x = myGUIdata.audioData;
fs = myGUIdata.samplingFrequency;
myGUIdata.audoPlayTimerPeriod = 0.057;
switch get(myGUIdata.player,'running')
    case 'off'
        myGUIdata.player = audioplayer(x/max(abs(x))*0.9,fs);
        set(myGUIdata.playOriginalButton,'string','STOP');
        guidata(handles.VTmainpulatorGUI,myGUIdata);
        set(myGUIdata.sgramCursorHandle,'visible','on');
        set(myGUIdata.player,'userdata',handles.VTmainpulatorGUI,'TimerFcn', ...
            @cursorUpdateTimerFunction,'TimerPeriod',myGUIdata.audoPlayTimerPeriod);
        set(myGUIdata.counterText,'visible','on');
        set(myGUIdata.lissajousPlotHandle,'visible','on');
        playblocking(myGUIdata.player);
        set(myGUIdata.playOriginalButton,'string','PLAY');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
    case 'on'
        stop(myGUIdata.player);
        set(myGUIdata.playOriginalButton,'string','PLAY');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
end;
end

function cursorUpdateTimerFunction(obj, event, string_arg)
handleForTimer = get(obj,'userData');
myGUIdata = guidata(handleForTimer);
currentPoint = get(myGUIdata.player,'CurrentSample');
fs = get(myGUIdata.player,'SampleRate');
index30ms = 1:round(0.03*fs);
TotalSamples = get(myGUIdata.player,'TotalSamples');
xx = myGUIdata.fundamentalComponent(min(TotalSamples,currentPoint+index30ms));
triggerPoint = find(xx(1:end-1).*xx(2:end)<0 & xx(1:end-1)<0, 1 );
ydata = myGUIdata.audioData(min(TotalSamples,currentPoint+index30ms+triggerPoint));
if sum(isnan(ydata))>0 || length(ydata) <1 || max(abs(ydata)) == 0
    set(myGUIdata.lissajousAxis,'ylim',[-1 1],'xlim',[1 index30ms(end)]);
else
    set(myGUIdata.lissajousAxis,'ylim',max(abs(ydata))*[-1 1],'xlim',[1 index30ms(end)]);
end;
set(myGUIdata.sgramCursorHandle,'xdata',currentPoint/fs*[1 1]);
set(myGUIdata.lissajousPlotHandle,'xdata',index30ms,'ydata',ydata);
set(myGUIdata.counterText,'string',[num2str(currentPoint/fs*1000,'%05.0f') '  ms']);
end

% --- Executes on button press in playSynthesizedButton.
function playSynthesizedButton_Callback(hObject, eventdata, handles)
% hObject    handle to playSynthesizedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
responseLengthInMs = 20;
fs = 8000;
lpcStructure = myGUIdata.lpcStructure;
f0Struct = myGUIdata.f0Struct;
%sgramStr = myGUIdata.sgramStr;
signalTime = (0:1/fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50);
for ii = 1:fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/fs*2*pi)).*(f0InSignalTime*ii<fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*fs):round(lpcStructure.frameShiftInMs/1000*fs))';
responseLength = round(responseLengthInMs/1000*fs);
fftLSynth = 2^ceil(log2(responseLength+length(baseIndex)+1));
olaIndex = (1:fftLSynth)';
for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStructure.preEmphasis*cos(theta(:))) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,1)*cos(theta(:)) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,2)*cos(2*theta(:));
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
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs)))+response;
end;
x = outputBuffer;
fs = 8000;
myGUIdata.synthesizedSignal = x;
myGUIdata.audoPlayTimerPeriod = 0.057;
set(myGUIdata.saveModifiedButton,'enable','off');
switch get(myGUIdata.player,'running')
    case 'off'
        myGUIdata.player = audioplayer(x/max(abs(x))*0.9,fs);
        set(myGUIdata.playSynthesizedButton,'string','STOP');
        guidata(handles.VTmainpulatorGUI,myGUIdata);
        set(myGUIdata.sgramCursorHandle,'visible','on');
        set(myGUIdata.player,'userdata',handles.VTmainpulatorGUI,'TimerFcn', ...
            @cursorUpdateTimerFunctionForSynth,'TimerPeriod',myGUIdata.audoPlayTimerPeriod);
        set(myGUIdata.counterText,'visible','on');
        set(myGUIdata.lissajousPlotHandle,'visible','on');
        playblocking(myGUIdata.player);
        set(myGUIdata.playSynthesizedButton,'string','Play resynthesis');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
    case 'on'
        stop(myGUIdata.player);
        set(myGUIdata.playSynthesizedButton,'string','Play resynthesis');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
end;
end

function cursorUpdateTimerFunctionForSynth(obj, event, string_arg)
handleForTimer = get(obj,'userData');
myGUIdata = guidata(handleForTimer);
currentPoint = get(myGUIdata.player,'CurrentSample');
fs = get(myGUIdata.player,'SampleRate');
index30ms = 1:round(0.03*fs);
TotalSamples = get(myGUIdata.player,'TotalSamples');
ydata = myGUIdata.synthesizedSignal(min(TotalSamples,currentPoint+index30ms));
if sum(isnan(ydata))>0 || length(ydata) <1 || max(abs(ydata)) == 0
    set(myGUIdata.lissajousAxis,'ylim',[-1 1],'xlim',[1 index30ms(end)]);
else
    set(myGUIdata.lissajousAxis,'ylim',max(abs(ydata))*[-1 1],'xlim',[1 index30ms(end)]);
end;
set(myGUIdata.sgramCursorHandle,'xdata',currentPoint/fs*[1 1]);
set(myGUIdata.lissajousPlotHandle,'xdata',index30ms,'ydata',ydata);
set(myGUIdata.counterText,'string',[num2str(currentPoint/fs*1000,'%05.0f') '  ms']);
end

% --- Executes on button press in saveModifiedButton.
function saveModifiedButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveModifiedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
outFileNameWav = ['modFile' datestr(now,30) '.wav'];
[file,path] = uiputfile(outFileNameWav,'Save modified sound and data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;
audiowrite([path file],myGUIdata.synthesizedSignal/max(abs(myGUIdata.synthesizedSignal))*0.9,8000);
outStatusFileName = [file(1:end-4) '.mat'];
modificationStructure.originalSignal = myGUIdata.audioData;
modificationStructure.samplingFrequency = myGUIdata.samplingFrequency;
modificationStructure.dataDescription = myGUIdata.dataDescription;
modificationStructure.lpcStructure = myGUIdata.lpcStructure;
modificationStructure.f0Struct = myGUIdata.f0Struct;
modificationStructure.sectionModifier = myGUIdata.sectionModifier;
modificationStructure.synthesizedSignal = myGUIdata.synthesizedSignal;
modificationStructure.samplingFrequencyOfSynth = 8000;
modificationStructure.f0ModificationRatio = myGUIdata.f0ModRatio;
modificationStructure.rawCompositeModifier = myGUIdata.rawCompositeModifier;
modificationStructure.rawSectionModifier = myGUIdata.rawSectionModifier;
modificationStructure.knobCoefficient = myGUIdata.knobCoefficient;
modificationStructure.shapeMagnifier = myGUIdata.shapeMagnifier;
modificationStructure.vtlMagnifier = myGUIdata.vtlMagnifier;
modificationStructure.whenSaved = datestr(now);
save([path outStatusFileName],'modificationStructure');
end

% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
close(myGUIdata.VTmainpulatorGUI);
end

% --- Executes on button press in snapShotButton.
% ---- modify July/17/2014 yoshimoto

function snapShotButton_Callback(hObject, eventdata, handles)
% hObject    handle to snapShotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
fileName = ['guiSnap' datestr(now,30) '.png'];

[file,path] = uiputfile(fileName,'Save the captured data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;

im = getframe(handles.VTmainpulatorGUI);
imwrite(im.cdata,[path file],'PNG');

end


% --- Executes on button press in playModifiedButton.
function playModifiedButton_Callback(hObject, eventdata, handles)
% hObject    handle to playModifiedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
responseLengthInMs = 20;
fs = 8000;
lpcStructure = myGUIdata.lpcStructure;
f0Struct = myGUIdata.f0Struct;
logArea = myGUIdata.lpcStructure.logAreaMatrix(1,:);
nSection = length(logArea);
locationList = (1:nSection)-1;
modifyerLocation = get(myGUIdata.compositePlotHandle,'xdata');
modifyerValue = get(myGUIdata.compositePlotHandle,'ydata');
sectionModifier = interp1(modifyerLocation*(nSection-1),modifyerValue,locationList,'linear','extrap');
myGUIdata.sectionModifier = sectionModifier;
%sgramStr = myGUIdata.sgramStr;
signalTime = (0:1/fs:f0Struct.temporalPositions(end))';
excitation = signalTime*0;
f0InSignalTime = interp1(f0Struct.temporalPositions,f0Struct.F0,signalTime,'linear',50)*myGUIdata.f0ModRatio;
for ii = 1:fs/40/2
    excitation = excitation+cos(cumsum(ii*f0InSignalTime/fs*2*pi)).*(f0InSignalTime*ii<fs/2);
end;
outputBuffer = excitation*0;
theta = lpcStructure.rawFAxis/fs*2*pi;
fftl = length(theta);
baseIndex = (-round(lpcStructure.frameShiftInMs/1000*fs):round(lpcStructure.frameShiftInMs/1000*fs))';
responseLength = round(responseLengthInMs/1000*fs);
fftLSynth = 2^ceil(log2(responseLength+length(baseIndex)+1));
olaIndex = (1:fftLSynth)';%myGUIdata.vtlMagnifier
frequencyAxis = (0:fftl-1)'/fftl*8000;
for ii = 1:length(lpcStructure.temporalPosition)
    inversePreprocessingShape = -2*log(1-lpcStructure.preEmphasis*cos(theta(:))) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,1)*cos(theta(:)) ...
        +myGUIdata.lpcStructure.cepstrumList(ii,2)*cos(2*theta(:));
    logArea = lpcStructure.logAreaMatrix(ii,:);
    newLogArea = logArea(:)+sectionModifier(:)/4;
    newRef = area2ref(exp(newLogArea));                                     %%newRef is reflection coefficient
    newalp = k2alp(newRef);                                                 %%newalp is predictor
    LPCspectrumLn = -2*log(abs(fft(newalp(:),fftl)));                       %%calculate LPC spectrum
    %LPCspectrumLn = -2*log(abs(fft(lpcStructure.lpcMatrix(ii,:)',fftl)));
    fixedLPCspectrum = LPCspectrumLn(:)+inversePreprocessingShape;
    fixedHalfLPCspectrum = interp1(frequencyAxis,fixedLPCspectrum, ...
        myGUIdata.vtlMagnifier*frequencyAxis(1:fftl/2));
    fixedLPCspectrum = [fixedHalfLPCspectrum;fixedHalfLPCspectrum(end-1:-1:2)];
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
    outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs))) = ...
        outputBuffer(min(length(excitation),olaIndex+round(lpcStructure.temporalPosition(ii)*fs)))+response;
end;
x = outputBuffer;
fs = 8000;
myGUIdata.synthesizedSignal = x;
myGUIdata.audoPlayTimerPeriod = 0.057;
set(myGUIdata.saveModifiedButton,'enable','on');
switch get(myGUIdata.player,'running')
    case 'off'
        myGUIdata.player = audioplayer(x/max(abs(x))*0.9,fs);
        set(myGUIdata.playModifiedButton,'string','STOP');
        guidata(handles.VTmainpulatorGUI,myGUIdata);
        set(myGUIdata.sgramCursorHandle,'visible','on');
        set(myGUIdata.player,'userdata',handles.VTmainpulatorGUI,'TimerFcn', ...
            @cursorUpdateTimerFunctionForSynth,'TimerPeriod',myGUIdata.audoPlayTimerPeriod);
        set(myGUIdata.counterText,'visible','on');
        set(myGUIdata.lissajousPlotHandle,'visible','on');
        playblocking(myGUIdata.player);
        set(myGUIdata.playModifiedButton,'string','Play modified');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
    case 'on'
        stop(myGUIdata.player);
        set(myGUIdata.playModifiedButton,'string','Play modified');
        set(myGUIdata.counterText,'visible','off');
        set(myGUIdata.lissajousPlotHandle,'visible','off');
        set(myGUIdata.sgramCursorHandle,'visible','off');
end;
end

% --- Executes on button press in saveOriginalButton.
function saveOriginalButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveOriginalButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
outFileName = ['rec' datestr(now,30) '.wav'];
[file,path] = uiputfile(outFileName,'Save the captured data');
if length(file) == 1 && length(path) == 1
    if file == 0 || path == 0
        %okInd = 0;
        disp('Save is cancelled!');
        return;
    end;
end;
audiowrite([path file],myGUIdata.audioData,myGUIdata.samplingFrequency);
end

% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
for focusIndex = 1:myGUIdata.numberOfBasisFunctions
    xdata = get(myGUIdata.controlHandleBundle(focusIndex).componentHandle,'xdata');
    ydata = xdata*0;
    set(myGUIdata.controlHandleBundle(focusIndex).componentControlHandle,'ydata',0);
    set(myGUIdata.controlHandleBundle(focusIndex).componentHandle,'ydata',ydata);
end;
myGUIdata.rawCompositeModifier = ydata(:);
set(myGUIdata.compositePlotHandle,'ydata',(myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
%set(myGUIdata.compositePlotHandle,'ydata',ydata);
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
myGUIdata = displayModifiedLPCspectrum(myGUIdata);
fftl = length(myGUIdata.lpcStructure.rawFAxis);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
guidata(handles.VTmainpulatorGUI,myGUIdata);
end

% --- Executes on button press in resetMagnifierButton.
function resetMagnifierButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetMagnifierButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
set(myGUIdata.shapeMagnifierHandle,'ydata',myGUIdata.shapeMagnifierReference);
myGUIdata.shapeMagnifier = 1;
set(myGUIdata.shapeMagnifierText,'string','1.00');
set(myGUIdata.shapeMagnifierFocusHandle,'ydata',myGUIdata.shapeMagnifierReference);
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);
set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifier*(myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier));
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
guidata(handles.VTmainpulatorGUI,myGUIdata);
end


% --- Executes on button press in resetVTLButton.
function resetVTLButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetVTLButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
set(myGUIdata.vtlHandle,'xdata',myGUIdata.vtlReference);
myGUIdata.vtlMagnifier = 1;
set(myGUIdata.vtlText,'string','1.00');
set(myGUIdata.vtlFocusHandle,'xdata',myGUIdata.vtlReference);
xdata = get(myGUIdata.originalVTAHandle,'xdata');
set(myGUIdata.manipulatedVTAHandle,'xdata',xdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
%set(myGUIdata.compositePlotHandle,'ydata',myGUIdata.shapeMagnifierReference*myGUIdata.rawCompositeModifier);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
guidata(handles.VTmainpulatorGUI,myGUIdata);
end


% --- Executes on button press in f0HandleResetButton.
function f0HandleResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to f0HandleResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
set(myGUIdata.F0Handle,'ydata',myGUIdata.F0ForDisplay);
set(myGUIdata.f0RatioText,'string','1.00');
myGUIdata.f0ModRatio = 1;
guidata(handles.VTmainpulatorGUI,myGUIdata);
end


% --- Executes on button press in resetKnobButton.
function resetKnobButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetKnobButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGUIdata = guidata(handles.VTmainpulatorGUI);
myGUIdata.knobCoefficient = myGUIdata.knobCoefficient*0;
myGUIdata = calculateSectionModifier(myGUIdata);
set(myGUIdata.compositePlotHandle,'ydata',...
    (myGUIdata.rawCompositeModifier+myGUIdata.rawSectionModifier)*myGUIdata.shapeMagnifier);
myGUIdata = updateKnobs(myGUIdata);
myGUIdata = displayAreaFunctions(myGUIdata,myGUIdata.xinit);
if get(myGUIdata.autoRegressiveModelButton,'value') == 1
    myGUIdata = displayModifiedLPCspectrum(myGUIdata);
    fftl = length(myGUIdata.lpcStructure.rawFAxis);
    set(myGUIdata.sgramImage,'xdata',myGUIdata.lpcStructure.temporalPosition,'ydata', ...
        myGUIdata.lpcStructure.rawFAxis(1:fftl/2+1)/1000, ...
        'cdata',max(-80,myGUIdata.fixedLPCspectrumInDB));
    set(myGUIdata.spectrogramAxis,'ylim',[0 4],'xlim',...
        [myGUIdata.lpcStructure.temporalPosition(1) myGUIdata.lpcStructure.temporalPosition(end)]);
end;
guidata(handles.VTmainpulatorGUI,myGUIdata);
end
