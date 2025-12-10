function rangemin_event_frames=detect_rangemin(range, metadata, localframe, filename)

% usage: rangemin_event_frames=detect_rangemin(range, metadata, localframe, filename);
% detects rangemin events using a peak-finding algorithm on inverted range
% rangemin events are closest-approach events, kind of like intercept or approach or contact

% helpful reminder about some of matlab's findpeaks parameters:
% Threshold - Minimum height difference between a peak and its neighbors
% MinPeakDistance — Minimum peak separation in time
% MinPeakWidth — Minimum peak width

%params:
MinPeakProminence=10;
MinPeakDistance=100;
MinPeakHeight=-10; %minimum range = 10
MinPeakWidth=32;
Threshold=0;

tic
fprintf('\ndetecting range minima... ')
rangemin_event_frames=[];
region=1:length(range); %you could make it shorter for param tuning

[pks_raw,rangemin_raw,w_raw,p_raw]  = findpeaks(-range(region), 'MinPeakProminence',MinPeakProminence, 'MinPeakDistance',MinPeakDistance, 'MinPeakHeight',MinPeakHeight, 'Threshold',Threshold, 'MinPeakWidth', MinPeakWidth);
%clean jumps that are just cricket drops, or flickers (p too high)
for f=rangemin_raw'
    cricketdrop = metadata{contains(metadata.filename, filename{f}), 'cricketdrop'};
    captureframe = metadata{contains(metadata.filename, filename{f}), 'captureframe'};
    if localframe(f)>cricketdrop+200 & localframe(f)<captureframe         %keep frames after cricket drop (+100 for travel) and before capture
            rangemin_event_frames=[rangemin_event_frames; f];
            % pks=[pks; pks_raw(find(f==rangemin_raw))];
            % w=[w; w_raw(find(f==rangemin_raw))];
            % p=[p; p_raw(find(f==rangemin_raw))];
    end
end
fprintf(' done (%.0f sec)', toc)
fprintf('\n%d raw rangemin events found', length(rangemin_raw))
fprintf('\n%d cleaned rangemin events found', length(rangemin_event_frames))
