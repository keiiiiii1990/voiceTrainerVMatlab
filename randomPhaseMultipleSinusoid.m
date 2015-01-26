function y = randomPhaseMultipleSinusoid(f0,fs,nHarmonics,duration)

tt = (0:1/fs:duration)';
y = tt*0;
for ii = 1:nHarmonics
    y = y+sin(2*ii*pi*f0*tt+2*pi*rand);
end;
y = sqrt(2)/sqrt(nHarmonics)*y;
return;