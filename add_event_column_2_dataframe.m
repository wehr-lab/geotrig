% add_event_column_2_dataframe          
% will write new dataframe, having added a last column 
% with the event in each row 
% modelled after working_script_kip01 (after working_script from Mike)
% can cooncatenate two input datafiles, but doesn't require them

close all

% conditions to be compared (only used here in constructing Kip's FNs)
condition1 = 10003;
condition2 = 10013;

%load geos from dataframe
dataroot= 'X:\PreyCapture\A1Suppression';
dataframeFN = strcat('geos_condition_', num2str(condition1) ,'p0_wfames.csv');
metadataFN = strcat('metadata_condition_', num2str(condition1) ,'p0');
newdataframeFN = strcat('EVENTS_',dataframeFN);

% optionally, enter a second dataframe (e.g. with and w/o laser) to be concatenated
dataframeFN2 = strcat('geos_condition_', num2str(condition2) ,'p0_wfames.csv');
metadataFN2 = strcat('metadata_condition_', num2str(condition2) ,'p0');

% where to store new dataframe:
outputrootdir= 'X:\PreyCapture\A1Suppression\geo-trig-analysis-output';

if ~exist('dataframeFN2','var')
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata, dataframe] = geotrig_load_dataframe_kip01(dataroot, dataframeFN, metadataFN);
else
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata, dataframe] = geotrig_load_dataframe_kip01(dataroot, dataframeFN, metadataFN, dataframeFN2, metadataFN2);
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

% make and populate new column - you could easily change these event acronyms
% note: this assumes there are no double labels!
fprintf('\nAdding new column called "event"\n')
dataframe.event(:)= "";  % initializes the column to have all empty ("") strings
dataframe.event(cricket_jump_event_frames) = 'cj';
dataframe.event(rangemin_event_frames) = 'rm';
dataframe.event(contact_gain_event_frames) = 'cg';
dataframe.event(contact_loss_event_frames) = 'cl';
dataframe.event(target_loss_event_frames) = 'tl';
dataframe.event(chase_start_frames) = 'cs';
dataframe.event(chase_end_frames) = 'ce';
dataframe.event(pause_start_frames) = 'ps';
dataframe.event(pause_end_frames) = 'pe';
dataframe.event(wander_start_frames) = 'ws';
dataframe.event(wander_end_frames) = 'we';
dataframe.event(stalk_start_frames) = 'ss';
dataframe.event(stalk_end_frames) = 'se';
dataframe.event(approach_start_frames) = 'as';
dataframe.event(approach_end_frames) = 'ae';

% save table as csv
fprintf('\nSaving new CSV to %s\n',fullfile(outputrootdir,newdataframeFN))
writetable(dataframe, fullfile(outputrootdir,newdataframeFN));
