function cricket_jump_event_frames=detect_cricketjump(cricket_spd, metadata, localframe, filename)

% usage:cricket_jump_event_frames=detect_cricketjump(cricket_spd, metadata, localframe, filename)
% detects cricketjump events using a peak-finding algorithm on cricket speed

% helpful reminder about some of matlab's findpeaks parameters:
% Threshold - Minimum height difference between a peak and its neighbors
% MinPeakDistance — Minimum peak separation in time
% MinPeakWidth — Minimum peak width

%params:
cricket_jum_win=10;
MinPeakProminence=10;
MinPeakDistance=10;
MinPeakHeight=30;
MinPeakWidth=10;
Threshold=0;
p_flicker_thresh=100; %peak prominence higher than this is probably a DLC glitch

fprintf('\ndetecting cricket jumps... ')
cricket_jump_event_frames=[];
pks=[];w=[];p=[];
region=1:length(cricket_spd); %you could make it shorter for param tuning
[pks_raw,cricket_jump_raw,w_raw,p_raw]  = findpeaks(cricket_spd(region), 'MinPeakProminence',MinPeakProminence, 'MinPeakDistance',MinPeakDistance, 'MinPeakHeight',MinPeakHeight, 'Threshold',Threshold, 'MinPeakWidth', MinPeakWidth);

%clean jumps that are just cricket drops, or DLC glitches (p too high)
for f=cricket_jump_raw'
    cricketdrop = metadata{contains(metadata.filename, filename{f}), 'cricketdrop'};
    captureframe = metadata{contains(metadata.filename, filename{f}), 'captureframe'};
    if localframe(f)>cricketdrop+100 & localframe(f)<captureframe         %keep frames after cricket drop (+100 for travel) and before capture
        if p_raw(find(f==cricket_jump_raw))<p_flicker_thresh
            cricket_jump_event_frames=[cricket_jump_event_frames; f];
            %pks=[pks; pks_raw(find(f==cricket_jump_raw))];
            %w=[w; w_raw(find(f==cricket_jump_raw))];
            %p=[p; p_raw(find(f==cricket_jump_raw))];
        end
    end
end
fprintf('\n%d raw cricket jump events found', length(cricket_jump_raw))
fprintf('\n%d cleaned cricket jump events found', length(cricket_jump_event_frames))

%note that pks, w, p are not used, but you could inspect them if you want to tweak the peakfinding params
% pks: height of each peak
% w: width of each peak
% p: prominence of each peak

% % inspect cricket jumps by uncommenting the following
% figure
% hold on
% t=1:length(cricket_spd(region));
% plot(t, cricket_spd(region), 'linewidth', 2)
% plot(cricket_jump_event_frames, MinPeakHeight+0*cricket_jump_event_frames, '*')
% zoom xon
% shg
% for i=1:10
%     xlim([cricket_jump_event_frames(i)-400 cricket_jump_event_frames(i)+400]);
%     pause(.5)
% end

