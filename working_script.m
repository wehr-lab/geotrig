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


% approach, intercept  ... not sure if I want to include jen's version? see /Users/wehr/Documents/Analysis/Prey-Capture/preycapture_simple.m


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


%get laser statistics for a single event type
[event_counts_ON, event_counts_OFF, pvalue, RateRatio_ON_vs_OFF] = GetEventLaserStats(rangemin_event_frames, 'rangemin', metadata, filename, 1);

%get laser statistics for lots of event types
fprintf('\nLaser statistics for events')
clear event_counts_ON event_counts_OFF
for e=1:length(event_types)
    event_frames=eval(event_types{e});
    [event_counts_ON(e), event_counts_OFF(e), pvalue(e), RateRatio_ON_vs_OFF(e)] = GetEventLaserStats(event_frames, event_types{e}, metadata, filename, 0);
    fprintf('\n%s: %d ON, %d OFF, ON ratio=%.1f, p=%.3f', event_types{e}, event_counts_ON(e), event_counts_OFF(e), RateRatio_ON_vs_OFF(e), pvalue(e))
end
% note since there are precisely the same number of contact gain/loss events, the stats for them are identical

%bar graph of event laser statistics computed above
figure
b= barh([event_counts_ON; event_counts_OFF]', 'grouped');
set(b(1), 'FaceColor', 'c')
set(b(2), 'FaceColor', 'k')
yticks(1:length(event_types))
yticklabels(event_types)
set(gca, 'TickLabelInterpreter', 'none') %prevents _ in labels from being interpreted as a LaTex-style subscript
for i=1:length(pvalue)
    if pvalue(i)<.05 %should probably do a multiple comparisons correction here
        text( max(event_counts_ON(i), event_counts_OFF(i))+100, i, '*', 'fontsize', 24)
    end
end
set(gcf, 'pos', [-1276 -1050  1154  2383])
set(gca, 'fontsize', 18)
xlabel('event count')
th=text(1, 1, '* p<0.05 laser effect on event rate, using Poisson regression', 'fontsize', 12)
set(th, 'units', 'normalized', 'position', [.6 .025])

%export stats to csv file using fprintf
fid=fopen(fullfile(outputrootdir, 'event_counts.csv'), 'w');
fprintf(fid, 'event_type, event_counts_ON, event_counts_ON, RateRatioONvsOFF, p-value');
for e=1:length(event_types)
    fprintf(fid, '\n%s, %d, %d, %.4f,  %.4f', event_types{e}, event_counts_ON(e), event_counts_OFF(e), RateRatio_ON_vs_OFF(e), pvalue(e));
end
fclose(fid);

%export stats to csv file by creating a table and then exporting it
%does exactly the same thing as the fprintf method above, but is more extensible
T = table(event_types(:), ...                % create table T from variables. Use (:) to enforce a column vector
    event_counts_ON(:), ...             
          event_counts_OFF(:), ...            
          RateRatio_ON_vs_OFF(:), ...   
          pvalue(:), ...                     
          'VariableNames', ...                % Assign column headers
          {'event_type', ...
           'event_counts_ON', ...
           'event_counts_OFF', ...
           'RateRatioONvsOFF', ...
           'p_value'});
outputFilePath = fullfile(outputrootdir, 'event_counts_table.csv');
writetable(T, outputFilePath); % Export the Table to CSV, using the writetable function, which handles file opening/closing and formatting automatically

%export event_frames data to a CSV that you could import into a dataframe
frame_data_cell = cell(length(event_types), 1); %  Initialize a cell array to hold the frame vectors
% Loop through the event names and pull the corresponding variable data
for i = 1:length(event_types)
    current_name = event_types{i};
    current_frames = eval(current_name);     % Retrieve the data from the workspace using eval
    frame_data_cell{i} = current_frames; % Store the entire frame vector into one cell of the cell array    
end
% Create the final table
T = table(event_names, frame_data_cell, 'VariableNames', {'EventType', 'EventFrameIndices'});
outputFilePath = fullfile(outputrootdir, 'event_frames_table.csv');
writetable(T, outputFilePath); % Export the Table to CSV, using the writetable function

%plot average geometries for a specific event type
max_num_clips=1*16;
GenerateEventVideoClips(rangemin_event_frames, 'rangemin', filename, localframe, 200, max_num_clips, outputrootdir)
GenerateEventVideoClips(pause_start_frames, 'pause', filename, localframe, pause_durs, max_num_clips, outputrootdir)













