function k = area2ref(s)
%   Area to reflection coefficients
%   s   : cross sectional area
%   k   : reflection coefficients

%   by Hideki Kawahara

n = length(s)-1;
k = zeros(n,1);
for ii=1:n
    k(ii) = (s(ii+1)-s(ii))/(s(ii+1)+s(ii));
end;
return;
