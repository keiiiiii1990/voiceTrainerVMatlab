function alp = k2alp(k)
%   Reflection coefficients to predictor
%   k   : reflection coefficient
%   alp : predictor

%   by Hideki Kawahara

n = length(k);
a = zeros(n,1);
b = zeros(n,1);
a(1) = k(1);
for ii=2:n
    for jj = 1:ii-1
        b(jj) = a(jj)-k(ii)*a(ii-jj);
    end;
    a = b;
    a(ii) = k(ii);
end;
alp = [1;-a];
return;
    