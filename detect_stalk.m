function [stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, cricket_spd, range, az) 

% usage: [stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, range, az) 
% detects stalk state,  a continuous binary state variable, defined by speed/range/az in a window
% also returns stalk_start_frames and stalk_end_frames as events, and stalk_durs 

% params:
stalk_speed_thresh=[3 15]; %[min max]  cm/s
stalk_range_thresh= [5 40]; % [min max],  cm
stalk_az_thresh=30; %min, in degrees
stalk_cspeed_thresh=1; % cricket speed min cm/s
stalk_winsize=1*200; %in frames (seconds*200fps)
min_stalk_dur=1*200; %in frames (seconds*200fps)


tic
fprintf('\ndetecting stalk... ')

medfilt_mouse_spd = medfilt1(mouse_spd, stalk_winsize, 'omitnan');
medfilt_az = medfilt1(az, stalk_winsize, 'omitnan');
medfilt_range = medfilt1(range, stalk_winsize, 'omitnan');
medfilt_cricket_spd = medfilt1(cricket_spd, stalk_winsize, 'omitnan');
%medfilt is fast (<<1s) and by doing it here, we can use a stalk-detection specific winsize 

stalk=zeros(size(mouse_spd)); %initialize to zero

% 1. Speed condition
spd_condition = (medfilt_mouse_spd > stalk_speed_thresh(1)) & ...
    (medfilt_mouse_spd < stalk_speed_thresh(2));

% 2. Range condition 
range_condition = (medfilt_range > stalk_range_thresh(1)) & ...
                  (medfilt_range < stalk_range_thresh(2));

% 3. Azimuth condition
az_condition = medfilt_az > stalk_az_thresh;

% 4. cricket_speed condition
cspeed_condition = medfilt_cricket_spd > stalk_cspeed_thresh;


% Combine all conditions (Boolean array)
stalk = spd_condition & range_condition & az_condition & cspeed_condition & cricket_present;

 fprintf(' done (%.1f sec)', toc)

stalk_start_frames=find(diff(stalk)==1);
stalk_end_frames=find(diff(stalk)==-1);
if stalk(end)==1
    stalk_end_frames=[stalk_end_frames f];
end
if length(stalk_start_frames) ~= length(stalk_end_frames) error('mismatched stalk start/stop'), end
stalk_durs=stalk_end_frames-stalk_start_frames;
fprintf('\nfound %d stalks, min duration %d frames (mean %.0f)', length(stalk_durs), min(stalk_durs),  mean(stalk_durs))

%exclude short stalks (optional)
keepidx=find(stalk_durs>=min_stalk_dur);
stalk_start_frames=stalk_start_frames(keepidx);
stalk_end_frames=stalk_end_frames(keepidx);
stalk_durs=stalk_durs(keepidx);
fprintf('\nafter excluding stalks <%d frames, kept %d stalks, min duration %d frames (mean %.0f)', min_stalk_dur, length(stalk_durs),  min(stalk_durs), mean(stalk_durs))
