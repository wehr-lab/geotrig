function [pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd, metadata, filename, localframe) 

% usage: [pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd, metadata, filename, localframe) 
% detects pause state,  a continuous binary state variable, defined by median-filtered speed below a threshold
% also returns pause_start_frames and pause_end_frames as events, and pause_durs 

% params:
pause_speed_thresh=2; %max,  cm/s
pause_winsize=.5*200; %in frames (seconds*200fps)
min_pause_dur=.5*200; %in frames (seconds*200fps)


tic
fprintf('\ndetecting pause... ')

medfilt_mouse_spd = medfilt1(mouse_spd, pause_winsize, 'omitnan');
%medfilt is fast (<<1s) and by doing it here, we can use a pause-detection specific winsize 

pause=zeros(size(mouse_spd)); %initialize to zero

% Speed condition 
spd_condition = medfilt_mouse_spd < pause_speed_thresh;

% exclude pauses at the end of the trial that are capture-related (optional)
exclude_frames=zeros(size(mouse_spd)); %initialize to zero
if (0)
    tic
    exclude_capture_winsize=10*200; %in frames (seconds*200fps)
    fprintf('\n excluding terminal capture-related pauses less than %d s before the final capture...', exclude_capture_winsize/200)
    for i=1:height(metadata)
        captureframe = metadata{i, 'captureframe'};
        fname=metadata{i, 'filename'};
        trialframes=find(contains(filename, fname));
        exclude_frames(trialframes)  = (localframe(trialframes) > captureframe-exclude_capture_winsize) & (localframe(trialframes) < captureframe)  ;
        %exclude frames that are less N s before captureframe and are on this trial
    end
    fprintf(' done (%.0f sec)\n', toc)
else
    fprintf('\n not excluding any terminal pauses')
end

% Combine  conditions 
pause = spd_condition & cricket_present & ~exclude_frames;

 fprintf(' done (%.1f sec)', toc)

pause_start_frames=find(diff(pause)==1);
pause_end_frames=find(diff(pause)==-1);
if pause(end)==1
    pause_end_frames=[pause_end_frames f];
end
if length(pause_start_frames) ~= length(pause_end_frames) error('mismatched pause start/stop'), end
pause_durs=pause_end_frames-pause_start_frames;
fprintf('\nfound %d pauses, min duration %d frames (mean %.0f)', length(pause_durs), min(pause_durs),  mean(pause_durs))

%exclude short pauses (optional)
keepidx=find(pause_durs>=min_pause_dur);
pause_start_frames=pause_start_frames(keepidx);
pause_end_frames=pause_end_frames(keepidx);
pause_durs=pause_durs(keepidx);
fprintf('\nafter excluding pauses <%d frames, kept %d pauses, min duration %d frames (mean %.0f)', min_pause_dur, length(pause_durs),  min(pause_durs), mean(pause_durs))




