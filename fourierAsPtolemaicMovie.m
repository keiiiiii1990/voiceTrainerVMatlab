function okInd = fourierAsPtolemaicMovie(sceneCycles,circleSizeMatrix,circlePhaseMatrix,makeMovie)
%   okInd = fourierAsPtolemaicMovie(sceneCycles,circleSizeMatrix,circlePhaseMatrix,makeMovie)
%
%   Example sctipt to show transition from sinusoid, square wave and
%   triangular wave
%
%sceneCycles = [2 2 2 2];
%circleSizeMatrix = [1 0 0 0 0; ...
%    1 0 1/3 0 0; ...
%    1 0 1/3 0 1/5; ...
%    1 0 1/3/3 0 1/5/5];
%circlePhaseMatrix = [0 0 0 0 0; ...
%    [1 0 1 0 1]*pi/2; ...
%    [1 0 1 0 1]*pi/2; ...
%    [1 0 1 0 1]*pi];
%
%makeMovie = 1; % 0: only display, 1: saves a series of pitcures in a directory
%okInd = fourierAsPtolemaicMovie(sceneCycles,circleSizeMatrix,circlePhaseMatrix,makeMovie);
%

%  Fourier analysis as Ptolemy
%   by Hideki Kawahara
%   07/Feb./2012
%   12/Feb./2012
%   This function was originally written by Hideki Kawahara.
%   This function is in public domain. Anything can be done with this.
%   The author does not have any responsibility of possible damages and
%   defficiencies, which this function may cause.

okInd = [];
figure;
th = 0:2*pi/60:2*pi;
yc = cos(-th);
xc = sin(-th);

%  get default figure size

figurePosition = get(gcf,'position');
circleHandle = subplot(211);
circleAxisPosition = get(circleHandle,'position');
subplot211Aspect = figurePosition.*circleAxisPosition;
subplot211Aspect = subplot211Aspect(4)/subplot211Aspect(3);

if makeMovie
directoryName = ['movie' datestr(now,30)];
mkdir(directoryName);
end;

%---- you can edit following three variables to make your own new movie
if 1 == 2
sceneCycles = [4 1 1 1 2 2 2 2 1 0 0 0 4];
circleSizeMatrix = [1 0 0 0 0; ...
    1 1 0 0 0; ...
    1 1 1 0 0; ...
    1 1 1 1 0; ...
    1 1 1 1 1; ...
    1 1 1 1 1; ...
    1 1 1 1 1; ...
    1 1 1 1 1; ...
    1 1 1 1 1; ...
    1 1 1 1 0; ...
    1 1 1 0 0; ...
    1 1 0 0 0; ...
    1 0 0 0 0];
circlePhaseMatrix = [0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    [1 1 1 1 1]*pi/2; ...
    0 0 0 0 0; ...
    [1 0 1 0 1]*pi/2; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0; ...
    0 0 0 0 0];
end;
np = length(th);
waveBuffer.shape = zeros(np*sum(sceneCycles)*2,1)*NaN;
waveBuffer.time = zeros(np*sum(sceneCycles)*2,1);
deltaT = 3/np;
waveBuffer.time(:,1) = ((1:size(waveBuffer.time))*deltaT)'*(-1);

kk = 1;
stableState = 1;
sceneLength = sceneCycles(1);
timeIndex = 1;
lastMaxAmplitude = sum(circleSizeMatrix(1,:));
while kk <= length(sceneCycles)
    switch stableState
        case 1
            sceneLength = sceneCycles(kk);
        case 0
            sceneLength = 1;
    end;
    for ii = 1:sceneLength*np
        switch stableState
            case 1
                circleSize = circleSizeMatrix(kk,:);
                circlePhase = circlePhaseMatrix(kk,:);
            case 0
                if kk == length(sceneCycles);break;end;
                circleSize = circleSizeMatrix(kk+1,:)*ii/np+ ...
                    circleSizeMatrix(kk,:)*(np-ii)/np;
                circlePhase = circlePhaseMatrix(kk+1,:)*ii/np+ ...
                    circlePhaseMatrix(kk,:)*(np-ii)/np;
        end;
        nHarmonics = length(circleSize);
        maxAmplitude = sum(abs(circleSize));
        subplot(circleHandle);
        xBase = 0; yBase = 0;
        for jj = 1:nHarmonics
            plot(xc*circleSize(jj)+xBase,yc*circleSize(jj)+yBase);
            hold on
            initialPhaseIndex = round(circlePhase(jj)/2/pi*np);
            if initialPhaseIndex<0;initialPhaseIndex = initialPhaseIndex+np;end;
            harmonicIndex = rem((ii-1)*jj+initialPhaseIndex,np)+1;
            nextxBase = xBase+xc(harmonicIndex)*circleSize(jj);
            nextyBase = yBase+yc(harmonicIndex)*circleSize(jj);
            plot([xBase nextxBase],[yBase nextyBase],'k','linewidth',2);
            xBase = nextxBase; yBase = nextyBase;
        end;
        axis(([0 2 -subplot211Aspect subplot211Aspect]*1.2)/subplot211Aspect*maxAmplitude+(-maxAmplitude)*[1 1 0 0])
        currentAxis = axis;
        waveBuffer.shape(timeIndex) = nextyBase;
        displayTime = (waveBuffer.time-waveBuffer.time(timeIndex))*maxAmplitude+maxAmplitude*1.1;
        maxDisplay = max(abs(waveBuffer.shape(displayTime<currentAxis(2))));
        if maxDisplay>lastMaxAmplitude
            maxAmplitude = maxDisplay;
        elseif lastMaxAmplitude>sum(abs(circleSize))
            maxAmplitude = lastMaxAmplitude-(lastMaxAmplitude-sum(abs(circleSize)))/(np/1.5);
        end;
        plot([maxAmplitude*1.1 100],0*[-100 100],'c','linewidth',4);
        hold on
        plot(0*[-100 100]+maxAmplitude*1.1,[-100 100],'c','linewidth',4);
        displayTime = (waveBuffer.time-waveBuffer.time(timeIndex))*maxAmplitude+maxAmplitude*1.1;
        plot(displayTime,waveBuffer.shape,'linewidth',2);
        plot([nextxBase maxAmplitude*1.1],[1 1]*nextyBase,'r','linewidth',2);
        axis(([0 2 -subplot211Aspect subplot211Aspect]*1.2)/subplot211Aspect*maxAmplitude+(-maxAmplitude)*[1 1 0 0]);
        axis off; hold off
        subplot(223);
        stem(circleSize,'linewidth',2);
        axis([0.8 length(circleSize)+0.2 -0.1 1.1]);
        xlabel('harmonic number');
        title('amplitude');
        subplot(224)
        stem(circlePhase,'linewidth',2);
        axis([0.8 length(circleSize)+0.2 [-1.1 1.1]*pi]);
        xlabel('harmonic number');
        title('phase');
        lastMaxAmplitude = maxAmplitude;
        drawnow
        if makeMovie
        print('-dpng','-r150',[directoryName '/frr' num2str(timeIndex,'%05d') '.png']);
        end;
        timeIndex = timeIndex+1;
    end;
    switch stableState
        case 1
            stableState = 0;
        case 0
            stableState = 1;
            kk = kk+1;
    end;
end;
okInd = 1;
return;