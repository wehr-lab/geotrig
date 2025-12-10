%example working script

% clear
close all

%load geos from dataframe
dataroot= '/Volumes/Projects/PreyCapture/ZIActivation';
dataframefilename = 'geos_allmice_alltrials_wfames.csv';
outputrootdir= '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output';
[az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(dataroot, dataframefilename);

cricket_present=get_cricket_present_frames(metadata, localframe, num_geoframes, filename);

%detect events
cricket_jump_event_frames=detect_cricketjump(cricket_spd, metadata, localframe, filename);
rangemin_event_frames=detect_rangemin(range, metadata, localframe, filename);
[contact, contact_gain_event_frames, contact_loss_event_frames]=detect_contact(range, cricket_present);
[target_loss_event_frames]=detect_target_loss(range, az, contact, metadata, localframe, num_geoframes, filename);
[chase, chase_start_frames, chase_end_frames, chase_durs]=detect_chase(cricket_present, mouse_spd, range, az);
[pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd);
[wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az);
[stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, cricket_spd, range, az);


% approach ? ... not sure if I want to include jen's version


%plot average geometries for a specific event type
closefigs=0;
delete(fullfile(outputrootdir, 'event-avg-geo.pdf')) %optional
plot_avg_geometries(cricket_jump_event_frames, 'cricket_jump_event_frames', mouse_spd, range, az, cricket_spd, closefigs, outputrootdir)

event_types={'cricket_jump_event_frames', 'rangemin_event_frames', 'contact_gain_event_frames', 'contact_loss_event_frames','target_loss_event_frames',...
    'chase_start_frames','pause_start_frames','pause_end_frames','wander_start_frames'}; %list of event types

%plot average geometries for a whole list of event types
for e=1:length(event_types)
    event_frames=eval(event_types{e});
    plot_avg_geometries(event_frames, event_types{e}, mouse_spd, range, az, cricket_spd, closefigs, outputrootdir)
end


%get laser statistics for events
fprintf('\nLaser statistics for events')
for e=1:length(event_types)
    event_frames=eval(event_types{e});
    [event_counts_ON, event_counts_OFF, pvalue, RateRatio_ON_vs_OFF] = GetEventLaserStats(event_frames, event_types{e}, metadata, filename, 0);
    fprintf('\n%s: %d ON, %d OFF, ON ratio=%.1f, p=%.3f', event_types{e}, event_counts_ON, event_counts_OFF, RateRatio_ON_vs_OFF, pvalue)
end
% note since there are precisely the same number of contact gain/loss events, the stats for them are identical

max_num_clips=1*16;
%plot average geometries for a specific event type
GenerateEventVideoClips(rangemin_event_frames, 'rangemin', filename, localframe, 200, max_num_clips, outputrootdir)
GenerateEventVideoClips(pause_start_frames, 'pause', filename, localframe, pause_durs, max_num_clips, outputrootdir)



















