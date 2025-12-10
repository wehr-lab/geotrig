function [contact, contact_gain_event_frames, contact_loss_event_frames]=detect_contact(range, cricket_present) 

% usage: [contact, contact_gain_event_frames, contact_loss_event_frames]=detect_contact(range, cricket_present); 
%
% detects contact, a continuous binary state variable, defined by range
% remaining below a threshold for a window of time
% contact_gain_event_frames are the frames when contact is gained 
% contact_loss_event_frames are the frames when contact is lost

% params:
winsize=.5*200; %in frames (seconds*200fps)
contact_thresh=10; %range definition of contact, cm


%% contact detection
tic
fprintf('\ndetecting contact... ')
contact=zeros(size(range));
medfilt_range = medfilt1(range, winsize, 'omitnan');

range_condition = medfilt_range < contact_thresh;
contact = range_condition & cricket_present;


contact_gain_event_frames=find(diff(contact)==1);
contact_loss_event_frames=find(diff(contact)==-1);

fprintf(' done (%.0f sec)', toc)

