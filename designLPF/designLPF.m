function [b] = designLPF(fp,fe,fs)
%   [b] = designLPF(fc,rp,fe,felevel,fs)
%   by Akagiri Hayato (2010.11.14)
%   (reference from "http://lis2.huie.hokudai.ac.jp/~toyo/MATLAB/#5-1")
%   
%   input
%   fp : passband edge frequency (Hz)
%   rp : passband ripple (dB)
%   fe : stopband edge frequency (Hz)
%   felevel : attenuation at stopband edge frequency (dB)
%   fs : sampling frequency (Hz)
%
%   output
%   b : coefficient at FIR filter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for LPF by MUSHRAM hidden anchor
%fp=3500;
%fe=4105;
%fs=44100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wp = fp/(fs/2);
we = fe/(fs/2);

width = we-wp;
N = ceil((6.6*pi)/(width*pi));
N = ceil((6.6*pi)/(width*pi))+mod(N,2);
N=N+4;
win=hamming(N);

wc = (wp+we)/2;
b=fir1(N-1,wc,win);

subplot(211);
plot(b);
legend('impulse responce');

[h,w]=freqz(b,1,512);
subplot(212);
plot((w/pi)*(fs/2),abs(h));
xlabel('Hz');
grid on;
legend('gain characteristic')

max_h=max(abs(h));
dB_h=20*log10(abs(h)/max_h);

figure;plot((dB_h));

rp=-min(dB_h(1:ceil(wp*512)+1))
felevel=-max(dB_h(ceil(we*512)+1:512))
feCenterlevel=-max(dB_h(ceil((4000/(fs/2))*512)+1:512))