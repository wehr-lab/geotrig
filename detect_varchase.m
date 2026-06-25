function [chase, chase_start_frames, chase_end_frames, chase_durs]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh, chase_range_thresh, chase_az_thresh) 

% usage: [chase, chase_start_frames, chase_end_frames, chase_durs]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh, chase_range_thresh, chase_az_thresh) 
% this is a modified version of detect_chase in which the thresholds are
% specified by the user rather than being fixed in the function definition.
% unlike detect_chase, each threshold is specified as [min max], you can use 0 or inf if you don't care about a particular endpoint
% the idea is to be able to detect a battery of different types of varchase
% (fast, medium, slow, near, medium, far, etc)
% detects varchase state,  a continuous binary state variable, defined by speed/range/az in a user specified window
% also returns chase_start_frames and chase_end_frames as events, and chase_durs 

% params:
chase_winsize=.25*200; %in frames (seconds*200fps)
min_chase_dur=.1*200; %in frames (seconds*200fps)


tic
fprintf('\ndetecting chase with mouse_spd in [%d %d] ... ', chase_speed_thresh(1), chase_speed_thresh(2))

 medfilt_mouse_spd = medfilt1(mouse_spd, chase_winsize, 'omitnan');
 medfilt_az = medfilt1(az, chase_winsize, 'omitnan');
 medfilt_range = medfilt1(range, chase_winsize, 'omitnan');

% medfilt_mouse_spd = medfilt1(mouse_spd, chase_winsize);
% medfilt_az = medfilt1(az, chase_winsize);
% medfilt_range = medfilt1(range, chase_winsize);

%medfilt is fast (<<1s) and by doing it here, we can use a chase-detection specific winsize 

chase=zeros(size(mouse_spd)); %initialize to zero

% 1. Speed condition
spd_condition = (medfilt_mouse_spd > chase_speed_thresh(1)) & ...
                (medfilt_mouse_spd < chase_speed_thresh(2));

% 2. Range condition
range_condition = (medfilt_range > chase_range_thresh(1)) & ...
                  (medfilt_range < chase_range_thresh(2));

% 3. Azimuth condition
az_condition = (medfilt_az > chase_az_thresh(1)) & ...
                (medfilt_az < chase_az_thresh(2));

% 4. Combine all conditions (Boolean array)
chase = spd_condition & range_condition & az_condition & cricket_present;

 fprintf(' done (%.1f sec)', toc)

chase_start_frames=find(diff(chase)==1);
chase_end_frames=find(diff(chase)==-1);
if chase(end)==1
    chase_end_frames=[chase_end_frames f];
end
if length(chase_start_frames) ~= length(chase_end_frames) error('mismatched chase start/stop'), end
chase_durs=chase_end_frames-chase_start_frames;
fprintf('\nfound %d chases, min duration %d frames (mean %.0f)', length(chase_durs), min(chase_durs),  mean(chase_durs))

%exclude short chases (optional)
keepidx=find(chase_durs>=min_chase_dur);
chase_start_frames=chase_start_frames(keepidx);
chase_end_frames=chase_end_frames(keepidx);
chase_durs=chase_durs(keepidx);
% regenerate chase from start/end frames, keeping only valid durations
chase=zeros(size(mouse_spd)); %initialize to zero
for i=1:length(chase_start_frames)
    chase(chase_start_frames(i):chase_end_frames(i))=1;
end

fprintf('\nafter excluding chases <%d frames, kept %d chases, min duration %d frames (mean %.0f)', min_chase_dur, length(chase_durs),  min(chase_durs), mean(chase_durs))
