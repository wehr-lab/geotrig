% plotTPM3
%quick duplication of PlotTPM but combining all search states into a single "search" composite
% and all pursuit states into a single "pursuit" to see
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
    [chase, chase_start_frames, chase_end_frames, chase_durs]=detect_chase(cricket_present, mouse_spd, range, az);

end


%% compile states into pursuit and search

pursuit_start_frames = [chase_start_frames; approach_start_frames];
%pursuit_start_frames = [chase_start_frames; approach_start_frames; contact_gain_event_frames; rangemin_event_frames];
search_start_frames = [pause_start_frames; wander_start_frames; stalk_start_frames];


%% make and plot TPMs

fig=figure;
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
        pursuit_start_frames(find(condition(pursuit_start_frames))), ...
        search_start_frames(find(condition(search_start_frames))), ...
        ...%contact_gain_event_frames(find(condition(contact_gain_event_frames))), ...
        cricket_jump_event_frames(find(condition(cricket_jump_event_frames))), ...
        target_loss_event_frames(find(condition(target_loss_event_frames))), ...
        rangemin_event_frames(find(condition(rangemin_event_frames))), ...
        };

    statenames={ ...
        'search' , ...
        'pursuit', ...
        ...%'contact', ...
        'cricketjump', ...
        'targetloss', ...
        'rangemin', ...
        };


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
    hTPM.Colormap=cool;
    hTPM.FontSize=14;

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
hTPM.FontSize=14;

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

%% plot tpm circle
opts.minProb=0.1;
opts.title='dark';
plot_tpm_circle(TPMdark, statenames, opts)
opts.title='light';
plot_tpm_circle(TPMlight, statenames, opts)


%% plot within-state shifts for the combined states
pursuit = [chase | approach];
search = [pause | wander | stalk ];

statenames={'pursuit', 'search'};
for statename=statenames

    state=eval(statename{:});
    dosqrt=1;

    figure
    tiledlayout(1,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 420])
    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
        histmatoff=sqrt(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    else
        histmaton=(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
        histmatoff=(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),(histmaton'-histmatoff'));
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    % x=getfield(thresholds, statename{:}, 'az');
    % y=getfield(thresholds, statename{:}, 'range');
    % pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    % rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    cl=clim;
    clim([-1 1]*max(abs(cl)))
    set(gca, 'fontsize', 18)

    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    else
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges),(histmaton-histmatoff));
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    % x=getfield(thresholds, statename{:}, 'range');
    % y=getfield(thresholds, statename{:}, 'mouse_spd');
    % pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    % rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    cl=clim;
    clim([-1 1]*max(abs(cl)))
set(gca, 'fontsize', 18)

    nexttile
    if dosqrt
        histmaton=(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
        histmatoff=(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    else
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),(histmaton-histmatoff));
    colormap blue_to_red
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    % x=getfield(thresholds, statename{:}, 'az');
    % y=getfield(thresholds, statename{:}, 'mouse_spd');
    % pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    % rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));
    cl=clim;
    clim([-1 1]*max(abs(cl)))
set(gca, 'fontsize', 18)
end