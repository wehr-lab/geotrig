function [intercept_event_frames]=detect_intercepts(approach_end_frames, range)

% usage: [intercept_event_frames]=detect_intercepts(approach_end_frames, range);
% detects intercept event, which is an approach that ends with range < 5cm
% (note this is the converse of failed approach, which is an approach that ends
% with range >5cm)
% note that you have to run detect_approach first to get approach_end_frames



% params:
range_thresh=5;    %max,  cm

intercept_event_frames=approach_end_frames(range(approach_end_frames)<=range_thresh);


fprintf('found %d intercepts\n', length(intercept_event_frames))
