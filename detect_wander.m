function [wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az) 

% usage: [wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az) 
% detects wander state,  a continuous binary state variable, defined by speed/range/az in a window
% also returns wander_start_frames and wander_end_frames as events, and wander_durs 

% params:
wander_speed_thresh=[5 15]; %[min max]  cm/s
wander_range_thresh= 15; % min,  cm
wander_az_thresh=30; %min, in degrees
wander_winsize=1*200; %in frames (seconds*200fps)
min_wander_dur=1*200; %in frames (seconds*200fps)


tic
fprintf('\ndetecting wander... ')

medfilt_mouse_spd = medfilt1(mouse_spd, wander_winsize, 'omitnan');
medfilt_az = medfilt1(az, wander_winsize, 'omitnan');
medfilt_range = medfilt1(range, wander_winsize, 'omitnan');
%medfilt is fast (<<1s) and by doing it here, we can use a wander-detection specific winsize 

wander=zeros(size(mouse_spd)); %initialize to zero

% 1. Speed condition
spd_condition = (medfilt_mouse_spd > wander_speed_thresh(1)) & ...
    (medfilt_mouse_spd < wander_speed_thresh(2));

% 2. Range condition 
range_condition = (medfilt_range > wander_range_thresh);
                  
% 3. Azimuth condition
az_condition = medfilt_az > wander_az_thresh;

% 4. Combine all conditions (Boolean array)
wander = spd_condition & range_condition & az_condition & cricket_present;

 fprintf(' done (%.1f sec)', toc)

wander_start_frames=find(diff(wander)==1);
wander_end_frames=find(diff(wander)==-1);
if wander(end)==1
    wander_end_frames=[wander_end_frames f];
end
if length(wander_start_frames) ~= length(wander_end_frames) error('mismatched wander start/stop'), end
wander_durs=wander_end_frames-wander_start_frames;
fprintf('\nfound %d wanders, min duration %d frames (mean %.0f)', length(wander_durs), min(wander_durs),  mean(wander_durs))

%exclude short wanders (optional)
keepidx=find(wander_durs>=min_wander_dur);
wander_start_frames=wander_start_frames(keepidx);
wander_end_frames=wander_end_frames(keepidx);
wander_durs=wander_durs(keepidx);
fprintf('\nafter excluding wanders <%d frames, kept %d wanders, min duration %d frames (mean %.0f)', min_wander_dur, length(wander_durs),  min(wander_durs), mean(wander_durs))
