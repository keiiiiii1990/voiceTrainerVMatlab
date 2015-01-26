function y = filterLPF(b,x,fs)

y = filter(b,1,x);

sound(y,fs);