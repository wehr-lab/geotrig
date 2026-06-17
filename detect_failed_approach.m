function [failed_approach_event_frames]=detect_failed_approach(approach_end_frames, range)

% usage: [failed_approach_event_frames]=detect_failed_approach(approach_end_frames, range);
% detects failed approach event, which is an approach that ends with range > 5cm
% (note this is the converse of intercept, which is an approach that ends
% with range <5cm)
% note that you have to run detect_approach first to get approach_end_frames



% params:
range_thresh=5;    %minimum,  cm

failed_approach_event_frames=approach_end_frames(range(approach_end_frames)>range_thresh);


fprintf('found %d failed approaches\n', length(failed_approach_event_frames))
