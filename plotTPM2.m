% plotTPM2
%quick duplication of PlotTPM but with only a single chase state, to see
%what it looks like
%

%% load geos from dataframe
if 0
    dataroot= '/Volumes/Projects/PreyCapture/ZIActivation';
    dataframefilename = 'geos_allmice_alltrials_wfames.csv';
    outputrootdir= '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output';
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(dataroot, dataframefilename);
    cd '/Users/wehr/Documents/Projects/Molly Shallow'

    cricket_present=get_cricket_present_frames(metadata, localframe, num_geoframes, filename);
    [lightdark light dark] = get_lightdark(filename, localframe, num_geoframes, metadata);
    [laseron] = get_laseron(filename, localframe, num_geoframes, metadata);

end

if 0

    %% detect events
    cricket_jump_event_frames=detect_cricketjump(cricket_spd, metadata, localframe, filename);
    rangemin_event_frames=detect_rangemin(range, metadata, localframe, filename);
    [contact, contact_gain_event_frames, contact_loss_event_frames]=detect_contact(range, cricket_present);
    [target_loss_event_frames]=detect_target_loss(range, az, contact, metadata, localframe, num_geoframes, filename);
    [pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd);
    [wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az);
    [stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, cricket_spd, range, az);
    [approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az);

end

%% detect battery of chase states

% %slow & far
% chase_speed_thresh1=[5 10]; %minimum,  cm/s
% chase_range_thresh1= [20 30]; % [min max] cm
% chase_az_thresh1=[30 45] ; %maximum, in degrees
%
% %medium
% chase_speed_thresh2=[10 30]; %minimum,  cm/s
% chase_range_thresh2= [10 20]; % [min max] cm
% chase_az_thresh2=[20 30] ; %maximum, in degrees
%
% %fast & and close
% chase_speed_thresh3=[30 100]; %minimum,  cm/s
% chase_range_thresh3= [0 10]; % [min max] cm
% chase_az_thresh3=[0 20] ; %maximum, in degrees
% Spacing using a power law (e.g., p=0.5 for square root spacing)
p = 1;
start_val=45;
end_val=0;
n=1;
azvalues = linspace(start_val^p, end_val^p, n+1) .^ (1/p);
p=1;
start_val=10;
end_val=40;
speedvalues = linspace(start_val^p, end_val^p, n+1) .^ (1/p);
p=1;
start_val=30;
end_val=0;
rangevalues = linspace(start_val^p, end_val^p, n+1) .^ (1/p);


for i=1:n
    %1 = slow, far, wide
    %n = fast, close, narrow

    % chase_speed_threshname=sprintf('chase_speed_thresh%d', i);
    % eval([chase_speed_threshname, '=[7*(i-1) 7*i];'])
    % chase_range_threshname=sprintf('chase_range_thresh%d', i);
    % eval([chase_range_threshname, '=[6*(10-i) 6*(11-i)];'])
    % chase_az_threshname=sprintf('chase_az_thresh%d', i);
    % eval([chase_az_threshname, '=[18*(10-i) 18*(11-i)];'])

    chase_speed_threshname=sprintf('chase_speed_thresh%d', i);
    eval([chase_speed_threshname, '=[speedvalues(i) speedvalues(i+1)];'])
    chase_range_threshname=sprintf('chase_range_thresh%d', i);
    eval([chase_range_threshname, '=[rangevalues(i+1) rangevalues(i)];'])
    chase_az_threshname=sprintf('chase_az_thresh%d', i);
    eval([chase_az_threshname, '=[azvalues(i+1) azvalues(i)];'])

    % %fix az = wide open
    % chase_az_threshname=sprintf('chase_az_thresh%d', i);
    % eval([chase_az_threshname, '=[azvalues(end) azvalues(1)];'])
    % %fix range = wide open
    % chase_range_threshname=sprintf('chase_range_thresh%d', i);
    % eval([chase_range_threshname, '=[rangevalues(end) rangevalues(1)];'])
    %fix range -> 2 windows
    % if i==1 %| i==2
    %     chase_range_threshname=sprintf('chase_range_thresh%d', i);
    %     eval([chase_range_threshname, '=[rangevalues(end-1) rangevalues(1)];'])
    % end

    fprintf('\nspeed %2d: %2.1f- %2.1f\t range: %2.1f- %2.1f\t az: %2.1f- %2.1f', i, eval(chase_speed_threshname), eval(chase_range_threshname), eval(chase_az_threshname))

end

for i=1:n
    str=sprintf('[chase%d, chase_start_frames%d, chase_end_frames%d, chase_durs%d]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh%d, chase_range_thresh%d, chase_az_thresh%d);', i, i, i, i, i, i, i);
    eval(str)
end

%
% [chase1, chase_start_frames1, chase_end_frames1, chase_durs1]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh1, chase_range_thresh1, chase_az_thresh1);
% [chase2, chase_start_frames2, chase_end_frames2, chase_durs2]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh2, chase_range_thresh2, chase_az_thresh2);
% [chase3, chase_start_frames3, chase_end_frames3, chase_durs3]=detect_varchase(cricket_present, mouse_spd, range, az, chase_speed_thresh3, chase_range_thresh3, chase_az_thresh3);


fig=figure
    tiledlayout(1,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1530 420])

%filter by condition
for c=1:2
   
    if c==1
        condition=dark; condition_name='dark';
    elseif c==2
        condition=light; condition_name='light';
    end
    % if c==1
    %     condition=~laseron.*dark; condition_name='laseroff.dark';
    % elseif c==2
    %     condition=laseron.*dark; condition_name='laseron.dark';
    % end

    % %all states
    %
    % states = {...
    %     approach_start_frames(find(condition(approach_start_frames))), ...
    %     chase_start_frames1(find(condition(chase_start_frames1))), ...
    %     chase_start_frames2(find(condition(chase_start_frames2))), ...
    %     chase_start_frames3(find(condition(chase_start_frames3))), ...
    %     chase_start_frames4(find(condition(chase_start_frames4))), ...
    %     rangemin_event_frames(find(condition(rangemin_event_frames))), ...
    %     contact_gain_event_frames(find(condition(contact_gain_event_frames))), ...
    %     cricket_jump_event_frames(find(condition(cricket_jump_event_frames))), ...
    %     stalk_start_frames(find(condition(stalk_start_frames))), ...
    %     pause_start_frames(find(condition(pause_start_frames))), ...
    %     target_loss_event_frames(find(condition(target_loss_event_frames))), ...
    %     wander_start_frames(find(condition(wander_start_frames)))};
    %
    % statenames={ ...
    %     'approach', ...
    %     'chase1' , ...
    %     'chase2' , ...
    %     'chase3' , ...
    %     'rangemin' , ...
    %     'contact' , ...
    %     'cricket jump' , ...
    %     'stalk' , ...
    %     'pause' , ...
    %     'target loss' , ...
    %     'wander' };

    % Subset of hand-picked States: Chase, approach, Pause, Wander, stalk
    states = {...
        ...%approach_start_frames(find(condition(approach_start_frames))), ...
        chase_start_frames1(find(condition(chase_start_frames1))), ...
        ...%chase_start_frames2(find(condition(chase_start_frames2))), ...
        ...%chase_start_frames3(find(condition(chase_start_frames3))), ...
        stalk_start_frames(find(condition(stalk_start_frames))), ...
        pause_start_frames(find(condition(pause_start_frames))), ...
        wander_start_frames(find(condition(wander_start_frames)))};

    statenames={ ...
        ...%'approach', ...
        'chase1' , ...
        ...%'chase2' , ...
        ...%'chase3' , ...
        'stalk' , ...
        'pause' , ...
        'wander' };


    % 1. Create a combined table of all events
    num_states = length(states); % Chase, Pause, Wander, etc.
    all_events = [];
    for i = 1:num_states % For each state
        starts = states{i};
        % Create array: [start_frame, state_ID]
        all_events = [all_events; starts, i*ones(length(starts), 1)];
    end

    % 2. Sort events by start frame to get the chronological sequence
    [~, idx] = sort(all_events(:, 1));
    state_sequence = all_events(idx, 2);

    %optional: eliminate Self-Transitions
    state_sequence(diff(state_sequence) ~= 0);

    T = zeros(num_states, num_states); % The raw count matrix

    % Count transitions
    for i = 1:length(state_sequence) - 1
        current_state = state_sequence(i);
        next_state = state_sequence(i+1);

        T(current_state, next_state) = T(current_state, next_state) + 1;
    end

    % Normalize to probabilities (rows sum to 1)
    TPM = T ./ sum(T, 2);

    % Absorbing States: If a state never transitions to another state, the sum
    % of that row will be zero, causing a division error. You can add a small
    % "pseudocount" (Laplace smoothing) to every cell in T (e.g., T = T + 0.01)
    % if you want to ensure the matrix is fully defined for sparse datasets.

    figure(fig)
     nexttile

    hTPM=heatmap(statenames, statenames, TPM, ...
        'Title', ['Behavioral Transition Probabilities (rows sum to 1), ', condition_name], ...
        'Colormap', parula, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
    xlabel('To State')
    ylabel('From State')

    figure;
    hT=heatmap(statenames, statenames, T, ...
        'Title', ['Behavioral Transition Counts, ', condition_name], ...
        'Colormap', parula, 'GridVisible', 'off');
    xlabel('To State')
    ylabel('From State')


    %hierarchical clustering

    % Calculate optimal leaf order for Rows (Y-axis)
    distRows = pdist(TPM, 'euclidean');
    linkRows = linkage(distRows, 'ward');
    orderRows = optimalleaforder(linkRows, distRows);

    % 3. Calculate optimal leaf order for Columns (X-axis)
    distCols = pdist(TPM', 'euclidean');
    linkCols = linkage(distCols, 'ward');
    orderCols = optimalleaforder(linkCols, distCols);

    % 4. Update the heatmap with new ordering
    % Convert numeric order to string cell arrays based on original labels
    hTPM.YDisplayData = hTPM.YData(orderRows);
    hTPM.XDisplayData = hTPM.XData(orderCols);
    title(hTPM, sprintf('Behavioral Transition Probabilities\n(rows sum to 1), asymmetric clustering, %s', condition_name))
    set(gca, 'fontsize', 14)
    hTPM.Colormap=hot;

    %clustergram is pretty cool, but the object handles are complex and it's
    %very quirky to customize, you have to find lots of hidden handles
    %cgo=clustergram(TPM)
    %set(cgo,'RowLabels',statenames,'ColumnLabels',statenames)
    %xlabel('To State') %hard to do
    %ylabel('From State')
    %fig = findall(0, 'Type', 'Figure', 'Tag', 'Clustergram');
    % The main heatmap axes have the tag 'HeatMapAxes'
    %ax = findall(fig, 'Type', 'Axes', 'Tag', 'HeatMapAxes');
    % 3. Modify the FontSize directly on the axes
    %set(ax, 'FontSize', 14);
     
    condition_names{c}=condition_name;

    if c==1 
         orderRowsdark=orderRows;
        orderColsdark=orderCols;
        TPMdark=TPM;
        hTPMdark=hTPM;
    elseif c==2 
        TPMlight=TPM;
         orderRowslight=orderRows;
        orderColslight=orderCols;
        hTPMlight=hTPM;
    end
end

%% plot dark-light diff

figure(fig);
nexttile
hTPM=heatmap(statenames, statenames, TPMdark-TPMlight, ...
    'Title', sprintf('Behavioral Transition Probabilities (rows sum to 1), %s - %s', condition_names{1}, condition_names{2} ), ...
    'Colormap', redbluecmap, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
xlabel('To State')
ylabel('From State')
cl=clim;
clim(max(abs(cl))*[-1 1])
hTPM.YDisplayData = hTPM.YData(orderRowsdark);
hTPM.XDisplayData = hTPM.XData(orderColsdark);

hTPMlight.YDisplayData = hTPM.YData(orderRowsdark);
hTPMlight.XDisplayData = hTPM.XData(orderColsdark);

clmax=max([clim(hTPMlight) clim(hTPMdark)]);
clim(hTPMdark, [0 clmax])
clim(hTPMlight, [0 clmax])

% fig=figure
%     tiledlayout(1,3, "TileSpacing","compact")
%     set(gcf, "Position", [440 880 1300 420])
% % nexttile
% hTPMdark.Parent = fig;
% hTPMdark.Layout.Tile = 1;
% hTPMlight.Parent = fig;
% hTPMlight.Layout.Tile = 2;
% hTPM.Parent = fig;
% hTPM.Layout.Tile = 3;
% hTPM.Position = []; 
% hTPM.ActivePositionProperty = 'position';
% shg

