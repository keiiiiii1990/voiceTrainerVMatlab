function cleanedStr = cleaningF0RawTrajectoryV2(r)

%%

contiguousThreshold = 2.5;
minimumSegment = 20; % ms

peakLocation = r.candidates.exLocationList(1,:);
totalFrames = length(peakLocation);
contiguousIndicator = [0 abs(diff(peakLocation))<contiguousThreshold] & ...
    [abs(diff(peakLocation))<contiguousThreshold 0] & (r.rawPeriodicity'>0.7);
indexTable = 1:length(peakLocation);
onsets = indexTable(diff([0 contiguousIndicator])>0);
offsets = indexTable(diff([0 contiguousIndicator])<0);

minSegments = min(length(onsets),length(offsets));
offsets = offsets(1:minSegments);
onsets = onsets(1:minSegments);

segmentList = zeros(minSegments,2);
segmentID = 0;
for ii = 1:minSegments
    if offsets(ii)-onsets(ii) > minimumSegment
        segmentID = segmentID+1;
        segmentList(segmentID,:) = [onsets(ii) offsets(ii)];
    end;
end;
segmentList = segmentList(1:segmentID,:);

%%

vuv = r.f0Raw*0;
for ii = 1:segmentID
    vuv(segmentList(ii,1):segmentList(ii,2)) = 1;
end;

%%

reliableF0 = r.f0Raw;
reliablePeriodicityIndex = r.rawPeriodicity;
latestLocation = r.candidates.exLocationList(1,segmentList(1,1));
for ii = segmentList(1,1)-1:-1:1
    nCandidate = sum(r.candidates.exLocationList(:,ii)>0);
    [minDistance,minIndex] = min(abs(r.candidates.exLocationList(1:nCandidate,ii)-latestLocation));
    reliableF0(ii) = r.candidates.f0Initial(minIndex,ii);
    reliablePeriodicityIndex(ii) = r.candidates.periodicityList(minIndex,ii);
    latestLocation = r.candidates.exLocationList(minIndex,ii);
end;
%%
latestLocation = r.candidates.exLocationList(1,segmentList(1,2));
for ii = 1:segmentID-1
    for jj = segmentList(ii,2)+1:ceil((segmentList(ii,2)+segmentList(ii+1,1))/2)
        nCandidate = sum(r.candidates.exLocationList(:,jj)>0);
        [minDistance,minIndex] = min(abs(r.candidates.exLocationList(1:nCandidate,jj)-latestLocation));
        reliableF0(jj) = r.candidates.f0Initial(minIndex,jj);
        reliablePeriodicityIndex(jj) = r.candidates.periodicityList(minIndex,jj);
        latestLocation = r.candidates.exLocationList(minIndex,jj);
    end;
    latestLocation = r.candidates.exLocationList(1,segmentList(ii+1,1));
    for jj = segmentList(ii+1,1)-1:-1:floor((segmentList(ii,2)+segmentList(ii+1,1))/2)
        nCandidate = sum(r.candidates.exLocationList(:,jj)>0);
        [minDistance,minIndex] = min(abs(r.candidates.exLocationList(1:nCandidate,jj)-latestLocation));
        reliableF0(jj) = r.candidates.f0Initial(minIndex,jj);
        reliablePeriodicityIndex(jj) = r.candidates.periodicityList(minIndex,jj);
        latestLocation = r.candidates.exLocationList(minIndex,jj);
    end;
    latestLocation = r.candidates.exLocationList(1,segmentList(ii+1,2));
end;
for jj = segmentList(segmentID,2)+1:totalFrames
    nCandidate = sum(r.candidates.exLocationList(:,jj)>0);
    [minDistance,minIndex] = min(abs(r.candidates.exLocationList(1:nCandidate,jj)-latestLocation));
    reliableF0(jj) = r.candidates.f0Initial(minIndex,jj);
    reliablePeriodicityIndex(jj) = r.candidates.periodicityList(minIndex,jj);
    latestLocation = r.candidates.exLocationList(minIndex,jj);
end;
cleanedStr = r;
cleanedStr.vuv = vuv;
cleanedStr.f0Raw = reliableF0(:);
cleanedStr.rawPeriodicity = reliablePeriodicityIndex;
cleanedStr.note = 'TrajectoryCleanedBycleaningF0RawTrajectory';
cleanedStr.cleanedAt = datestr(now);

