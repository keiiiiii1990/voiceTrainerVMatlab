%%  F0 tracking test
%   23/Sep./2012
%   by Hideki Kawahara

fs = 44100;
f0List = 32*2.^(0:1/6:log2(1025/32))';
duration = 10;
tt = (0:1/fs:duration)';
snrList = [20,30,40,50,60]';
f0RmsErrorList = zeros(length(f0List),length(snrList));
f0RmsErrorListRev = zeros(length(f0List),length(snrList));

for jj = 1:length(snrList)
    noiseAmplitude = 10.0^(-snrList(jj)/20);
    for ii = 1:length(f0List)
        f0 = f0List(ii);
        x = tt*0;
        for kk = 1:floor(fs/f0/2)
            x = x+sin(2*pi*f0*kk*tt+2*pi*rand);
        end;
        x = x+randn(length(x),1)*std(x)*noiseAmplitude;
        r = extractInitialF0bySymmetryICSP2013(x,fs);
        f0Raw = r.f0Raw(200:end-200);
        vuv = abs((f0Raw-f0)/f0)<0.2;
        refR = refineExGxOutputRev2Sqx(x,fs,r);
        f0Refined = refR.f0Refined(200:end-200);
        vuvRev = abs((f0Refined-f0)/f0)<0.2;
        f0RmsErrorListRev(ii,jj) = sqrt(mean((f0Refined(vuv)-f0).^2))/f0;
        f0RmsErrorList(ii,jj) = sqrt(mean((f0Raw(vuv)-f0).^2))/f0;
    end;
end;

%%

figure;loglog(f0List,f0RmsErrorList);grid on;
set(gca,'fontsize',15)
axis([30 1024 10^(-5) 10^(-1)]);
xlabel('frequency (Hz)');
ylabel('relative rms error');
title('initial estimate: flat sig. flat noise');
print -depsc initialF0rmsErrorFreq.eps

figure;loglog(f0List,f0RmsErrorListRev);grid on;
set(gca,'fontsize',15)
axis([30 1024 10^(-5) 10^(-1)]);
xlabel('frequency (Hz)');
ylabel('relative rms error');
title('refined estimate: flat sig. flat noise');
print -depsc refinedF0rmsErrorFreq.eps

%%  test with phase randomization

fs = 44100;
f0List = 32*2.^(0:1/6:log2(1025/32))';
duration = 1;
numberOfIterations = 10;
tt = (0:1/fs:duration)';
snrList = [20,30,40,50,60]';
f0RmsErrorList = zeros(length(f0List),length(snrList));
f0RmsErrorListRev = zeros(length(f0List),length(snrList));

for jj = 1:length(snrList)
    noiseAmplitude = 10.0^(-snrList(jj)/20);
    for ii = 1:length(f0List)
        f0 = f0List(ii);
        dataCount = 0;
        f0Raw = zeros(round(1000*duration*numberOfIterations),1);
        f0Refined = zeros(round(1000*duration*numberOfIterations),1);
        for ll = 1:numberOfIterations
            x = tt*0;
            for kk = 1:floor(fs/f0/2)
                x = x+sin(2*pi*f0*kk*tt+2*pi*rand);
            end;
            x = x+randn(length(x),1)*std(x)*noiseAmplitude;
            r = extractInitialF0bySymmetryICSP2013(x,fs);
            tmpf0Raw = r.f0Raw(100:end-100);
            vuv = abs((tmpf0Raw-f0)/f0)<0.2;
            refR = refineExGxOutputRev2Sqx(x,fs,r);
            tmpf0Refined = refR.f0Refined(100:end-100);
            f0Raw(dataCount+(1:sum(double(vuv)))) = tmpf0Raw(vuv);
            f0Refined(dataCount+(1:sum(double(vuv)))) = tmpf0Refined(vuv);
            dataCount = dataCount+sum(double(vuv));
        end;
        f0Raw = f0Raw(1:dataCount);
        f0Refined = f0Refined(1:dataCount);
        f0RmsErrorListRev(ii,jj) = sqrt(mean((f0Refined-f0).^2))/f0;
        f0RmsErrorList(ii,jj) = sqrt(mean((f0Raw-f0).^2))/f0;
    end;
end;

%%

figure;loglog(f0List,f0RmsErrorList);grid on;
set(gca,'fontsize',15)
axis([30 1024 10^(-5) 10^(-1)]);
xlabel('frequency (Hz)');
ylabel('relative rms error');
title('initial estimate: flat sig. flat noise');
legend('20dB','30dB','40dB','50dB','60dB','location','southeast');
print -depsc initialF0rmsErrorFreqRnd.eps

figure;loglog(f0List,f0RmsErrorListRev);grid on;
set(gca,'fontsize',15)
axis([30 1024 10^(-5) 10^(-1)]);
xlabel('frequency (Hz)');
ylabel('relative rms error');
title('refined estimate: flat sig. flat noise');
legend('20dB','30dB','40dB','50dB','60dB','location','northeast');
print -depsc refinedF0rmsErrorFreqRnd.eps
