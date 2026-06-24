% plotTPM

if 0
    %% load geos from dataframe
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
    [pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd, metadata, filename, localframe);

    [wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az);
    [stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, cricket_spd, range, az);
    [approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az);
    failed_approach_event_frames=detect_failed_approach(approach_end_frames, range);
    intercept_event_frames=detect_intercepts(approach_end_frames, range);
end


%% detect chase states by passing thresholds to detect_varchase

azvalues=[0 30];
rangevalues = [5 20];
speedvalues_hotpursuit = [30 50];
speedvalues_chase = [20 30];
speedvalues_follow = [5 20];

[hotpursuit, hotpursuit_start_frames, hotpursuit_end_frames, hotpursuit_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_hotpursuit, rangevalues, azvalues);
[chase, chase_start_frames, chase_end_frames, chase_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_chase, rangevalues, azvalues);
[follow, follow_start_frames, follow_end_frames, follow_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_follow, rangevalues, azvalues);




%% enforce hotpursuit/chase/follow wins for all overlaps with stalk/wander
hotpursuit_overlaps = hotpursuit & (stalk | wander);
chase_overlaps = chase & (stalk | wander) ;
follow_overlaps = follow & (stalk | wander) ;
fprintf('\nfound %d overlap frames between hotpursuit and stalk/wander (%.0f%% of hotpursuit frames)', sum(hotpursuit_overlaps), 100*sum(hotpursuit_overlaps)/sum(hotpursuit))
fprintf('\nfound %d overlap frames between chase and stalk/wander (%.0f%% of chase frames)', sum(chase_overlaps), 100*sum(chase_overlaps)/sum(chase))
fprintf('\nfound %d overlap frames between follow and stalk/wander (%.0f%% of follow frames)', sum(follow_overlaps), 100*sum(follow_overlaps)/sum(follow))
stalk(hotpursuit_overlaps)=0;
stalk(chase_overlaps)=0;
stalk(follow_overlaps)=0;
wander(hotpursuit_overlaps)=0;
wander(chase_overlaps)=0;
wander(follow_overlaps)=0;
fprintf('\nremoved all overlapping hotpursuit/chase/follow and wander/stalk frames, they are now exclusively hotpursuit/chase/follow')

% enforce stalk wins for overlaps with wander
stalk_wander_overlaps = wander & stalk ;
fprintf('\nfound %d overlap frames between wander and stalk (%.0f%% of wander frames, %.0f%% of stalk frames)', sum(stalk_wander_overlaps), 100*sum(stalk_wander_overlaps)/sum(wander), 100*sum(stalk_wander_overlaps)/sum(stalk))
wander(stalk_wander_overlaps)=0;
fprintf('\nremoved all overlapping wander frames, they are now exclusively stalk')

% enforce pause wins for all overlaps
pause_overlaps = pause & (stalk | wander) ;
fprintf('\nfound %d overlap frames between pause and wander OR stalk (%.0f%% of pause frames, %.0f%% of wander frames, %.0f%% of stalk frames)', sum(pause_overlaps), 100*sum(pause_overlaps)/sum(pause),  100*sum(pause_overlaps)/sum(wander), 100*sum(pause_overlaps)/sum(stalk))
wander(pause_overlaps)=0;
stalk(pause_overlaps)=0;
fprintf('\nremoved all overlapping stalk/wander frames, they are now exclusively pause')

pause_overlaps = pause & (hotpursuit | chase | follow ) ;
fprintf('\nfound %d overlap frames between pause and hotpursuit/chase/follow (%.0f%% of pause frames)', sum(pause_overlaps), 100*sum(pause_overlaps)/sum(pause))
hotpursuit(pause_overlaps)=0;
chase(pause_overlaps)=0;
follow(pause_overlaps)=0;
fprintf('\nremoved all overlapping hotpursuit/chase/follow frames, they are now exclusively pause')

%detect none of the above states
[noneoftheabove, noneoftheabove_start_frames, noneoftheabove_end_frames]=detect_noneoftheabove(hotpursuit, chase, follow, pause, stalk, wander);


if 0
%% save events and states to csv

% Each eventframes list → its own CSV column (pad with NaN)

    maxLen = max(cellfun(@numel, {...
        cricket_jump_event_frames, ...
        rangemin_event_frames, ...
        contact_gain_event_frames, contact_loss_event_frames, ...
        pause_start_frames, pause_end_frames, pause_durs, ...
        wander_start_frames, wander_end_frames, ...
        stalk_start_frames, stalk_end_frames, ...
        approach_start_frames, approach_end_frames, ...
        failed_approach_event_frames, ...
        intercept_event_frames, ...
        hotpursuit_start_frames, hotpursuit_end_frames, ...
        chase_start_frames, chase_end_frames, ...
        follow_start_frames, follow_end_frames, ...
        noneoftheabove_start_frames, noneoftheabove_end_frames}));
    pad = @(x) [x(:); NaN(maxLen - numel(x), 1)];
    eventsTable = table(...
        pad(cricket_jump_event_frames), ...
        pad(rangemin_event_frames), ...
        pad(contact_gain_event_frames),  ...
        pad(contact_loss_event_frames), ...
        pad(pause_start_frames),  ...
        pad(pause_end_frames), ...
        pad(wander_start_frames),  ...
        pad(wander_end_frames), ...
        pad(stalk_start_frames),  ...
        pad(stalk_end_frames), ...
        pad(approach_start_frames),  ...
        pad(approach_end_frames), ...
        pad(failed_approach_event_frames), ...
        pad(intercept_event_frames), ...
        pad(hotpursuit_start_frames),  ...
        pad(hotpursuit_end_frames), ...
        pad(chase_start_frames),  ...
        pad(chase_end_frames), ...
        pad(follow_start_frames),  ...
        pad(follow_end_frames), ...
        pad(noneoftheabove_start_frames) , ...
        pad(noneoftheabove_end_frames), ...
        'VariableNames', {
        'cricket_jump_event_frames', ...
        'rangemin_event_frames', ...
        'contact_gain_event_frames',  ...
        'contact_loss_event_frames', ...
        'pause_start_frames',  ...
        'pause_end_frames', ...
        'wander_start_frames',  ...
        'wander_end_frames', ...
        'stalk_start_frames',  ...
        'stalk_end_frames', ...
        'approach_start_frames',  ...
        'approach_end_frames', ...
        'failed_approach_event_frames', ...
        'intercept_event_frames', ...
        'hotpursuit_start_frames',  ...
        'hotpursuit_end_frames', ...
        'chase_start_frames',  ...
        'chase_end_frames', ...
        'follow_start_frames',  ...
        'follow_end_frames', ...
        'noneoftheabove_start_frames' , ...
        'noneoftheabove_end_frames'  });
    writetable(eventsTable, '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output/mike''s states and events/events.csv' );


    % State logicals → one rectangular CSV
    stateTable=table(contact(:),...
        pause(:),...
        wander(:),...
        stalk(:),...
        chase(:),...
        follow(:),...
        hotpursuit(:),...
        approach(:),...
        noneoftheabove(:), ...
        'VariableNames', {
        'contact',...
        'pause',...
        'wander',...
        'stalk',...
        'chase',...
        'follow',...
        'hotpursuit',...
        'approach',...
        'noneoftheabove', ...
}    );
    writetable(stateTable,    '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output/mike''s states and events/states.csv');
end


%% compute and plot TPM

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
   

   

    states = {...
        hotpursuit_start_frames(find(condition(hotpursuit_start_frames))), ...
        chase_start_frames(find(condition(chase_start_frames))), ...
        follow_start_frames(find(condition(follow_start_frames))), ...
        stalk_start_frames(find(condition(stalk_start_frames))), ...
        wander_start_frames(find(condition(wander_start_frames))), ...
        pause_start_frames(find(condition(pause_start_frames))), ...
        noneoftheabove_start_frames(find(condition(noneoftheabove_start_frames)))};
  
    num_states = length(states); % Chase, Pause, Wander, etc.

    statenames={ ...
        'hot pursuit' , ...
        'chase' , ...
        'following' , ...
        'stalk' , ...
        'wander' , ...
        'pause' , ...
        'none'};


    % 1. Create a combined table of all events
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
    % state_sequence = state_sequence(diff([0; state_sequence]) ~= 0);

    T = zeros(num_states, num_states); % The raw count matrix

    noneIdx=find(contains(statenames,'none'));

    % Count transitions
    for i = 1:length(state_sequence) - 1
        current_state = state_sequence(i);
        next_state = state_sequence(i+1);
        if current_state == noneIdx || next_state == noneIdx %skip counting none transitions
            continue;
        end
        T(current_state, next_state) = T(current_state, next_state) + 1;
    end

    % Absorbing States: If a state never transitions to another state, the sum
    % of that row will be zero, causing a division error. You can add a small
    % "pseudocount" (Laplace smoothing) to every cell in T (e.g., T = T + 0.01)
    % if you want to ensure the matrix is fully defined for sparse datasets.
    % T = T + 0.01;

    %exclude noneoftheabove state from further consideration
    keepidx=setdiff(1:num_states, noneIdx);
    T=T(keepidx, keepidx);
    statenames=statenames(keepidx);
    states=states(keepidx);
    num_states = length(states);


    % Normalize to probabilities (rows sum to 1)
    TPM = T ./ sum(T, 2);


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

    % bar graph of counts of each state
    clear statecounts
    for i = 1:num_states % For each state
        starts = states{i};
        statecounts(i) = length(starts);
    end
    figure
    bar(statenames, statecounts)
    ylabel('state counts')
    title(condition_name)


    %save results for plotting differences later, outside the loop
    condition_names{c}=condition_name;
    if c==1
        orderRowsdark=orderRows;
        orderColsdark=orderCols;
        TPMdark=TPM;
        hTPMdark=hTPM;
        statecountsdark=statecounts;
    elseif c==2
        TPMlight=TPM;
        orderRowslight=orderRows;
        orderColslight=orderCols;
        hTPMlight=hTPM;
        statecountslight=statecounts;
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



%% plot tpm circle
opts.minProb=0.1;
opts.colorbar=0;
opts.title='dark';
plot_tpm_circle(TPMdark, statenames, opts)
opts.title='light';
plot_tpm_circle(TPMlight, statenames, opts)

% plot_tpm_diff_circle(TPMdark-TPMlight, statenames, 'MinAbsDiff', 0.001);
% this still looks very messy

%% bar graph of dark-light state counts
figure
diffs = statecountsdark - statecountslight;
b = bar(statenames, diffs);
b.FaceColor = 'flat';
b.CData = (diffs' <= 0) * [0 0 1] + (diffs' > 0) * [1 0 0];
ylabel('dark-light state count diffs')
title('dark-light')

figure
b = bar(statenames, [ statecountsdark; statecountslight]);
ylabel(' state counts ')
title('dark vs light')
legend('dark', 'light', 'location', 'northwest')

%% plot eigenspectrum
plot_tpm_spectrum(TPMdark)

%% plot State probability trajectories from each starting state

%plot for a specific initial condition, such as pause
pi0=contains(statenames, 'pause'); %0s for all states except 1 at the match
numsteps=200;
plot_occupancy_over_time(TPMdark, numsteps,    'stateNames', statenames, 'pi0', pi0);

%plot for each  initial condition
for i=1:length(statenames)
    pi0=zeros(size(statenames));
    pi0(i)=1;
    plot_occupancy_over_time(TPMdark, numsteps,    'stateNames', statenames, 'pi0', pi0);
end

