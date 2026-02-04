%example working script (from Mike)
% working_script_kip01
% clear
close all

% conditions to be compared
condition1 = 10003;
condition2 = 10013;

%load geos from dataframe
dataroot= 'X:\PreyCapture\A1Suppression';
dataframeFN = strcat('geos_condition_', num2str(condition1) ,'p0_wfames.csv');
metadataFN = strcat('metadata_condition_', num2str(condition1) ,'p0');

% optionally, enter a second dataframe (e.g. with and w/o laser) to be concatenated
dataframeFN2 = strcat('geos_condition_', num2str(condition2) ,'p0_wfames.csv');
metadataFN2 = strcat('metadata_condition_', num2str(condition2) ,'p0');

outputrootdir= 'X:\PreyCapture\A1Suppression\geo-trig-analysis-output';

if ~exist('dataframeFN2')
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe_kip01(dataroot, dataframeFN, metadataFN);
else
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe_kip01(dataroot, dataframeFN, metadataFN, dataframeFN2, metadataFN2);
end

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
[approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az);

if 0
% approach, intercept  ... not sure if I want to include jen's version? see /Users/wehr/Documents/Analysis/Prey-Capture/preycapture_simple.m
    %%% define approaches
    approach = abs(az)<30 & mouse_spd>5;
    approach = medfilt1(approach,31); %%% removes brief periods and connects across gaps, on order of 0.5sec

    approachStarts = find(diff(approach)>0)+1;
    firstApproach = min(approachStarts); %first time point of approach
end

%plot average geometries for a specific event type
closefigs=0;
delete(fullfile(outputrootdir, 'event-avg-geo.pdf')) %optional
plot_avg_geometries(cricket_jump_event_frames, 'cricket_jump_event_frames', mouse_spd, range, az, cricket_spd, closefigs, outputrootdir)

event_types={'cricket_jump_event_frames', 'rangemin_event_frames', 'contact_gain_event_frames', 'contact_loss_event_frames','target_loss_event_frames',...
   'chase_start_frames','pause_start_frames','pause_end_frames','wander_start_frames'}; %list of event types

event_typesShort={'cricket_jump_event', 'rangemin_event', 'contact_gain_event', 'contact_loss_event','target_loss_event',...
    'chase_start','pause_start','pause_end','wander_start'}; %list of event types

%plot average geometries for a whole list of event types
for e=1:length(event_types)
    event_frames=eval(event_types{e});
    plot_avg_geometries(event_frames, event_types{e}, mouse_spd, range, az, cricket_spd, closefigs, outputrootdir)
end


%get statistics for a single event type
% totalnumframesON/OFF added fpr normalization
verbose=1;
[event_counts_ON, event_counts_OFF, pvalue, RateRatio_ON_vs_OFF, totalnumframesON, totalnumframesOFF] = GetEventLaserStats_kip02(rangemin_event_frames, 'rangemin', metadata, filename, condition1, condition2, verbose);

%get statistics for lots of event types
fprintf('statistics for events\n')
clear event_counts_ON event_counts_OFF
for e=1:length(event_types)
    event_frames=eval(event_types{e});
    [event_counts_ON(e), event_counts_OFF(e), pvalue(e), RateRatio_ON_vs_OFF(e), totalnumframesON, totalnumframesOFF] = GetEventLaserStats_kip02(event_frames, event_types{e}, metadata, filename, condition1, condition2, verbose);
    fprintf('%s: %d ON, %d OFF, ON ratio=%.1f, p=%.3f\n', event_types{e}, event_counts_ON(e), event_counts_OFF(e), RateRatio_ON_vs_OFF(e), pvalue(e))
end
% note since there are precisely the same number of contact gain/loss events, the stats for them are identical

%bar graph of event laser statistics computed above
% replaced event_counts_ON/OFF with normedEventsON/OFF
figure
normedEventsON = event_counts_ON/totalnumframesON*200;
normedEventsOFF = event_counts_OFF/totalnumframesOFF*200;

b= barh([normedEventsON; normedEventsOFF]', 'grouped');
Xlim = xlim;
Ylim = ylim;
set(b(1), 'FaceColor', 'c')
set(b(2), 'FaceColor', 'k')
yticks(1:length(event_types))
yticklabels(event_typesShort)
text(Xlim(2)+.01,Ylim(2)+.5,'n','fontsize',14)
text(Xlim(2)/2,Ylim(2)*.75,['nframes=' num2str(totalnumframesOFF)],'fontsize',12)
text(Xlim(2)/2,Ylim(2)*.7,['nframes=' num2str(totalnumframesON)],'color','blue','fontsize',12)

set(gca, 'TickLabelInterpreter', 'none') %prevents _ in labels from being interpreted as a LaTex-style subscript
for i=1:length(pvalue)
    text(Xlim(2)+.01, i-.25, num2str(event_counts_ON(i)), 'Color', 'blue', 'fontsize', 10)
    text(Xlim(2)+.01, i+.25, num2str(event_counts_OFF(i)), 'fontsize', 10)
    if pvalue(i)<.05 %should probably do a multiple comparisons correction here
        text( max(normedEventsON(i), normedEventsOFF(i))+.01, i, '*', 'fontsize', 24)
    end
end
%set(gcf, 'pos', [-1276 -1050  1154  2383])
set(gca, 'fontsize', 18)
xlabel('event count/sec')
th=text(1, 1, '* p<0.05', 'fontsize', 12)       %condition effect on event rate, using Poisson regression
set(th, 'units', 'normalized', 'position', [.6 -.025])

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
%this creates a table with N rows where N = number of event types, column 1
%is each event name, and columns 2-inf are all the event frames for that event type 
frame_data_cell = cell(length(event_types), 1); %  Initialize a cell array to hold the frame vectors
% Loop through the event names and pull the corresponding variable data
for i = 1:length(event_types)
    current_name = event_types{i};
    current_frames = eval(current_name);     % Retrieve the data from the workspace using eval
    frame_data_cell{i} = current_frames; % Store the entire frame vector into one cell of the cell array    
end
% Create the final table
T = table(event_types, frame_data_cell, 'VariableNames', {'EventType', 'EventFrameIndices'});
%T = table(event_names, frame_data_cell, 'VariableNames', {'EventType', 'EventFrameIndices'});
outputFilePath = fullfile(outputrootdir, 'event_frames_table.csv');
writetable(T, outputFilePath); % Export the Table to CSV, using the writetable function

%plot average geometries for a specific event type
max_num_clips=1*16;
GenerateEventVideoClips(rangemin_event_frames, 'rangemin', filename, localframe, 200, max_num_clips, outputrootdir)
GenerateEventVideoClips(pause_start_frames, 'pause', filename, localframe, pause_durs, max_num_clips, outputrootdir)













