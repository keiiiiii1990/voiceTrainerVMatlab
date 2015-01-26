function candidateStr = candidatesFromTwoMaps(r)
%temporalPositions = r.temporalPositions;
%numberOfLocations = length(temporalPositions);
y = r.updatedPeriodicityMap;
numberOfLocations = length(y);

maxFixedP = 5;
f0CandidatesMap = r.f0CandidatesMap;
periodicityList = zeros(maxFixedP,numberOfLocations);
exLocationList = zeros(maxFixedP,numberOfLocations);
f0Initial = zeros(maxFixedP,numberOfLocations);
xid = (1:length(r.fcList))';
dd = diff(r.updatedPeriodicityMap,1);
zx = zeros(1,size(dd,2));
uu = ([zx;dd].*[dd;zx]<0 & [dd;zx]<0);
for ii = 1:numberOfLocations
    if sum(uu(:,ii))==0
        [periodicityList(1,ii),exLocationList(1,ii)] = ...
            max(r.updatedPeriodicityMap(:,ii));
        f0Initial(1,ii) = f0CandidatesMap(exLocationList(1,ii),ii);
    else 
        uuc = uu(:,ii);
        a = (y(xid(uuc)-1,ii)+y(xid(uuc)+1,ii)-2*y(xid(uuc),ii))/2;
        b = (y(xid(uuc)+1,ii)-y(xid(uuc)-1,ii))/2;
        xFrag = -b./(a*2);
        x = xid(uuc)+xFrag;
        yx = y(xid(uuc),ii)+a.*xFrag.^2+b.*xFrag;
        [yxs,idx] = sort(yx,'descend');
        xs = x(idx);
        exLocationList(1:min(maxFixedP,length(x)),ii) = xs(1:min(maxFixedP,length(x)));
        periodicityList(1:min(maxFixedP,length(x)),ii) = yxs(1:min(maxFixedP,length(x)));
        f0c = f0CandidatesMap(xid(uuc),ii);
        f0m = f0CandidatesMap(xid(uuc)-1,ii);
        f0p = f0CandidatesMap(xid(uuc)+1,ii);
        f0i = (1-abs(xFrag)).*f0c+abs(xFrag).*f0m.*(xFrag<0)+abs(xFrag).*f0p.*(xFrag>0);
        f0is = f0i(idx);
        f0Initial(1:min(maxFixedP,length(x)),ii) = f0is(1:min(maxFixedP,length(x)));
    end;
end;
candidateStr.exLocationList = exLocationList;
candidateStr.periodicityList = periodicityList;
candidateStr.f0Initial = f0Initial;
return;