function [approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az)

% usage: [approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az)
% detects approach state,  a continuous binary state variable, defined by speed/az
% also returns approach_start_frames and approach_end_frames as events,
% approach_durs, and first_approach_frames (the first frame of the first approach for each trial)
% this uses jen's algorithm from Hoy et al 2016, as implemented in preycapture_simple.m

% from preycapture_simple.m (ran separately on each trial): 
% %% define approaches 
% approach = abs(az)<30 & spd>5;
% approach = medfilt1(approach,31); %%% removes brief periods and connects across gaps, on order of 0.5sec
%
% approachStarts = find(diff(approach)>0)+1;
% firstApproach = min(approachStarts); %first time point of approach


% params:
approach_speed_thresh=5; %minimum,  cm/s
approach_az_thresh=30; %maximum, in degrees
approach_winsize=.5*200; %in frames (seconds*200fps)
min_approach_dur=.25*200 %in frames (seconds*200fps)

tic
fprintf('\ndetecting approach... ')

approach=zeros(size(mouse_spd)); %initialize to zero

% Speed condition
spd_condition = mouse_spd > approach_speed_thresh;

% Azimuth condition
az_condition = az < approach_az_thresh;

% Combine all conditions (Boolean array)
approach = spd_condition & az_condition & cricket_present;

approach = medfilt1(single(approach),approach_winsize); %%% removes brief periods and connects across gaps, on order of 0.5sec

fprintf(' done (%.1f sec)', toc)

approach_start_frames=find(diff(approach)>0);
approach_end_frames=find(diff(approach)<0);
if approach(end)==1
    approach_end_frames=[approach_end_frames f];
end
if length(approach_start_frames) ~= length(approach_end_frames) error('mismatched approach start/stop'), end
approach_durs=approach_end_frames-approach_start_frames;
fprintf('\nfound %d approachs, min duration %d frames (mean %.0f)', length(approach_durs), min(approach_durs),  mean(approach_durs))

% get first approach frame for each trial, using cricket_present which indicates cricketdrop for each trial

first_approach_frames=[];
cricket_drop_frames=find(diff(cricket_present)==1);
for f=cricket_drop_frames(:)'
    c=    find(approach_start_frames>f, 1);
    first_approach_frames=[first_approach_frames approach_start_frames(c)];
end

%exclude short approachs (optional)
keepidx=find(approach_durs>=min_approach_dur);
approach_start_frames=approach_start_frames(keepidx);
approach_end_frames=approach_end_frames(keepidx);
approach_durs=approach_durs(keepidx);
fprintf('\nafter excluding approachs <%d frames, kept %d approachs, min duration %d frames (mean %.0f)', min_approach_dur, length(approach_durs),  min(approach_durs), mean(approach_durs))
