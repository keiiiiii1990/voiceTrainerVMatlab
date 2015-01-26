function s = ref2area(k)

n = length(k);
s = zeros(n+1,1);
s(end) = 1;
for ii = n:-1:1
    s(ii) = s(ii+1)*(1-k(ii))/(1+k(ii));
end;
