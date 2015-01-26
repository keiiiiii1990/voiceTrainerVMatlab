function refinedStr = refineExGxOutputRev2Sqx(x,fs,periodicityStructure,optStr)
%   refinedStr = refineExGxOutputRev2Sqx(x,fs,periodicityStructure)

%   by Hideki Kawahara
%   18/Oct./2011
%   23/Oct./2011 revised for option
%   13/Nov./2011 revised for HR tech meeting (x is test)
%   18/Aug./2012 safe guard version

numberOfHarmonics = 6;
windowStretch = 3;
numberOfInterference = 1;
if nargin == 4
    if isfield(optStr,'numberOfHarmonics')
        numberOfHarmonics = optStr.numberOfHarmonics;end;
    if isfield(optStr,'windowStretch')
        windowStretch = optStr.windowStretch;end;
    if isfield(optStr,'numberOfInterference')
        numberOfInterference = optStr.numberOfInterference;end;
end;
tic
refinedF0 = periodicityStructure.f0Raw;
for ii = 1:length(periodicityStructure.temporalPositions)
    rc = tandemGxF0RefineSqxRev(x,fs,periodicityStructure.temporalPositions(ii),...
        periodicityStructure.f0Raw(ii),numberOfHarmonics,windowStretch,numberOfInterference);
    refinedF0(ii) = rc.f0;
end;
optStr.numberOfHarmonics = numberOfHarmonics;
optStr.windowStretch = windowStretch;
optStr.numberOfInterference = numberOfInterference;
refinedStr = periodicityStructure;
refinedStr.f0Refined = refinedF0;
refinedStr.elapsedTimeForRefine = toc;
refinedStr.refineConditions = optStr;
return;