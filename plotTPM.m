% plotTPM

close all
clear results
if 0
    %% load geos from dataframe
    if ispc
        dataroot= '\\wehr-nas.uoregon.edu\Projects\PreyCapture\ZIActivation';
        dataframefilename = 'geos_allmice_alltrials_wfames.csv';
        outputrootdir= '\\wehr-nas.uoregon.edu\Projects\PreyCapture\ZIActivation\geo-trig-analysis-output';
        [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(dataroot, dataframefilename);
        cd 'C:\Users\wehrlab\Desktop\test_tpm'

    else
        dataroot= '/Volumes/Projects/PreyCapture/ZIActivation';
        dataframefilename = 'geos_allmice_alltrials_wfames.csv';
        outputrootdir= '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output';
        [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(dataroot, dataframefilename);
        cd '/Users/wehr/Documents/Projects/Molly Shallow'
    end

    cricket_present=get_cricket_present_frames(metadata, localframe, num_geoframes, filename);
    [lightdark light dark] = get_lightdark(filename, localframe, num_geoframes, metadata);
    [laseron] = get_laseron(filename, localframe, num_geoframes, metadata);

end

if 1

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
    %capture_frames=detect_capture(metadata, filename);
    %trial_durations=detect_trialdurs(metadata, filename);

end

fps=200;

%% detect chase states by passing thresholds to detect_varchase

azvalues=[0 30];
rangevalues_hotpursuit = [0 20];
rangevalues_chase = [5 20];
rangevalues_follow = [5 20];
speedvalues_hotpursuit = [25 100];
speedvalues_chase = [15 25];
speedvalues_follow = [5 15];

[hotpursuit, hotpursuit_start_frames, hotpursuit_end_frames, hotpursuit_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_hotpursuit, rangevalues_hotpursuit, azvalues);
[chase, chase_start_frames, chase_end_frames, chase_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_chase, rangevalues_chase, azvalues);
[follow, follow_start_frames, follow_end_frames, follow_durs]=detect_varchase(...
    cricket_present, mouse_spd, range, az, speedvalues_follow, rangevalues_follow, azvalues);




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

stalk_wander_overlaps = wander & stalk ;
fprintf('\nfound %d overlap frames between wander and stalk (%.0f%% of wander frames, %.0f%% of stalk frames)', sum(stalk_wander_overlaps), 100*sum(stalk_wander_overlaps)/sum(wander), 100*sum(stalk_wander_overlaps)/sum(stalk))

stalk_wins=1;
if stalk_wins
    % enforce stalk wins for overlaps with wander
    wander(stalk_wander_overlaps)=0;
    fprintf('\nstalk wins over wander: removed all overlapping wander frames, they are now exclusively stalk')
else
    % enforce wander wins for overlaps with stalk
    stalk(stalk_wander_overlaps)=0;
    fprintf('\nwander wins over stalk: removed all overlapping stalk/wander frames, they are now exclusively wander')
end

% enforce pause wins for all overlaps
pause_overlaps = pause & (stalk | wander) ;
fprintf('\nfound %d overlap frames between pause and wander OR stalk (%.0f%% of pause frames, %.0f%% of wander frames, %.0f%% of stalk frames)', sum(pause_overlaps), 100*sum(pause_overlaps)/sum(pause),  100*sum(pause_overlaps)/sum(wander), 100*sum(pause_overlaps)/sum(stalk))
wander(pause_overlaps)=0;
stalk(pause_overlaps)=0;
fprintf('\npause wins over stalk/wander: removed all overlapping stalk/wander frames, they are now exclusively pause')

pause_overlaps = pause & (hotpursuit | chase | follow ) ;
fprintf('\nfound %d overlap frames between pause and hotpursuit/chase/follow (%.0f%% of pause frames)', sum(pause_overlaps), 100*sum(pause_overlaps)/sum(pause))
hotpursuit(pause_overlaps)=0;
chase(pause_overlaps)=0;
follow(pause_overlaps)=0;
fprintf('\nremoved all overlapping hotpursuit/chase/follow frames, they are now exclusively pause')

%detect none of the above states
[noneoftheabove, noneoftheabove_start_frames, noneoftheabove_end_frames]=detect_noneoftheabove(hotpursuit, chase, follow, pause, stalk, wander);


% regenerate start and end frames now that overlaps are removed
hotpursuit_start_frames=find(diff(hotpursuit)==1);
chase_start_frames=find(diff(chase)==1);
follow_start_frames=find(diff(follow)==1);
stalk_start_frames=find(diff(stalk)==1);
wander_start_frames=find(diff(wander)==1);
pause_start_frames=find(diff(pause)==1);
noneoftheabove_start_frames=find(diff(noneoftheabove)==1);

hotpursuit_end_frames=find(diff(hotpursuit)==-1);
chase_end_frames=find(diff(chase)==-1);
follow_end_frames=find(diff(follow)==-1);
stalk_end_frames=find(diff(stalk)==-1);
wander_end_frames=find(diff(wander)==-1);
pause_end_frames=find(diff(pause)==-1);
noneoftheabove_end_frames=find(diff(noneoftheabove)==-1);
fprintf('\nregenerated start and end frames');
fprintf('\ntotal counts:')
fprintf('\nhotpursuit: %d', length( hotpursuit_end_frames));
fprintf('\nchase: %d', length( chase_end_frames));
fprintf('\nfollow: %d', length( follow_end_frames));
fprintf('\nstalk: %d', length( stalk_end_frames));
fprintf('\nwander: %d', length( wander_end_frames));
fprintf('\npause: %d', length( pause_end_frames));
fprintf('\nnoneoftheabove: %d', length( noneoftheabove_end_frames));








%% compute and plot TPM

fig=figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [440 440 1020 840])
figUC=figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [440 440 1020 840])

%filter by 2 conditions (dark & light)
% for c=1:2
% 
%     if c==1
%         condition=dark; condition_name='dark';
%     elseif c==2
%         condition=light; condition_name='light';
%     end

%filter by all 4 conditions (dark & light & laser on/off)
for c=1:4

    if c==1
        condition=laseron & dark; condition_name='dark laser on';
    elseif c==2
        condition=laseron & light; condition_name='light laser on';
    elseif c==3
        condition=~laseron & dark; condition_name='dark laser off';
    elseif c==4
        condition=~laseron & light; condition_name='light laser off';
    end

   

    states = {...
        hotpursuit_start_frames(find(condition(hotpursuit_start_frames))), ...
        chase_start_frames(find(condition(chase_start_frames))), ...
        follow_start_frames(find(condition(follow_start_frames))), ...
        stalk_start_frames(find(condition(stalk_start_frames))), ...
        wander_start_frames(find(condition(wander_start_frames))), ...
        pause_start_frames(find(condition(pause_start_frames))), ...
        ...%capture_frames(find(condition(capture_frames))), ...
        noneoftheabove_start_frames(find(condition(noneoftheabove_start_frames)))};
  
    num_states = length(states); % Chase, Pause, Wander, etc.

    statenames={ ...
        'hot pursuit' , ...
        'chase' , ...
        'following' , ...
        'stalk' , ...
        'wander' , ...
        'pause' , ...
        ...%'capture', ...
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
    captureIdx=find(contains(statenames,'capture'));

    % Count transitions
    for i = 1:length(state_sequence) - 1
        current_state = state_sequence(i);
        next_state = state_sequence(i+1);
        if current_state == noneIdx || next_state == noneIdx %skip counting none transitions
            continue;
        end
        if current_state == captureIdx  %skip counting transitions from capture to anything
            %(because capture terminates a trial)
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

        figure(figUC)
    nexttile

    hTPMuc=heatmap(statenames, statenames, TPM, ...
        'Title', ['Behavioral Transition Probabilities (unclustered), ', condition_name], ...
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
    bar((statenames), statecounts)
    ylabel('state counts')
    title(condition_name)


    %save results for plotting differences later, outside the loop
    condition_names{c}=condition_name;
    results(c).orderRows=orderRows;
    results(c).orderCols=orderCols;
    results(c).TPM=TPM;
    results(c).hTPM=hTPM;
    results(c).statecounts=statecounts;
    results(c).condition_name=condition_name

end

%% plot dark-light-laser diffs

figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
% 1-3 : dark on-off
% 2-4 : light on-off
% 1-2 : dark-light off
% 3-4 : dark-light on
A=[1 2 1 3];
B=[3 4 2 4];
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    hTPM=heatmap(statenames, statenames, results(a).TPM-results(b).TPM, ...
        'Title', sprintf('Difference TPM, %s - %s', results(a).condition_name, results(b).condition_name ), ...
        'Colormap', blue_to_red, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
    xlabel('To State')
    ylabel('From State')
    cl=clim;
    clim(max(abs(cl))*[-1 1])
    % cluster if desired:
    % hTPM.YDisplayData = hTPM.YData(results(3).orderRows);
    % hTPM.XDisplayData = hTPM.XData(results(3).orderCols);
    H(c)=hTPM;
end

%% plot tpm-circle for dark-light-laser diffs

figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
% 1-3 : dark on-off
% 2-4 : light on-off
% 1-2 : dark-light off
% 3-4 : dark-light on
A=[1 2 1 3];
B=[3 4 2 4];
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    diffMat=results(a).TPM-results(b).TPM;
    mytitle=sprintf('Difference TPM, %s - %s', results(a).condition_name, results(b).condition_name );
    plot_tpm_diff_circle(diffMat, statenames, 'Title', mytitle)

end

%%  repeat plot dark-light-laser diffs but with light-dark flipped

figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
% 1-3 : dark on-off
% 2-4 : light on-off
% 1-2 : dark-light off
% 3-4 : dark-light on
A=[1 2 2 4];
B=[3 4 1 3];
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    hTPM=heatmap(statenames, statenames, results(a).TPM-results(b).TPM, ...
        'Title', sprintf('Difference TPM, %s - %s', results(a).condition_name, results(b).condition_name ), ...
        'Colormap', blue_to_red, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
    xlabel('To State')
    ylabel('From State')
    cl=clim;
    clim(max(abs(cl))*[-1 1])
    % cluster if desired:
    % hTPM.YDisplayData = hTPM.YData(results(3).orderRows);
    % hTPM.XDisplayData = hTPM.XData(results(3).orderCols);
    H(c)=hTPM;
end

%% repeat plot  tpm-circle for dark-light-laser diffs but with light-dark flipped

figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
% 1-3 : dark on-off
% 2-4 : light on-off
% 1-2 : dark-light off
% 3-4 : dark-light on
A=[1 2 2 4];
B=[3 4 1 3];
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    diffMat=results(a).TPM-results(b).TPM;
    mytitle=sprintf('Difference TPM, %s - %s', results(a).condition_name, results(b).condition_name );
    plot_tpm_diff_circle(diffMat, statenames, 'Title', mytitle)

end





%plot without clustering
figure
tiledlayout(1,2, "TileSpacing","compact")
set(gcf, "Position", [440 440 1100 420])
nexttile
c=3; %dark laseroff
hTPMdark_uc=heatmap(statenames, statenames, results(c).TPM, ...
    'Title', sprintf('TPM %s, no clustering', results(c).condition_name), ...
    'Colormap', parula, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
xlabel('To State')
ylabel('From State')

nexttile
c=4; %light laseroff
hTPMdark_uc=heatmap(statenames, statenames, results(c).TPM, ...
    'Title', sprintf('TPM %s, no clustering', results(c).condition_name), ...
    'Colormap', parula, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
xlabel('To State')
ylabel('From State')
colormap hot;

%set all clustered heatmaps to the dark-off order
for c=1:4
    results(c).hTPM.YDisplayData=results(c).hTPM.YData(results(3).orderRows);
    results(c).hTPM.XDisplayData=results(c).hTPM.XData(results(3).orderCols);
end

if 0
    %% plot TPM heatmaps
    % (if you want to do it outside the loop)
    figure;
    tiledlayout(2,2, "TileSpacing","compact")
    set(gcf, "Position", [460 460 1020 840])
    for c=1:4
        nexttile
        h=heatmap(statenames, statenames, results(c).TPM, ...
            'Title', ['TPM, ', results(c).condition_name], ...
            'Colormap', hot, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
        xlabel('To State')
        ylabel('From State')

        h.YDisplayData=h.YData(results(3).orderCols);
        h.XDisplayData=h.XData(results(3).orderCols);

    end
end

%% plot tpm circle
opts.minProb=0.1;
opts.colorbar=0;
for c=1:4
    opts.title=results(c).condition_name;
    plot_tpm_circle(results(c).TPM, statenames, opts)
end

% opts.minProb=0.1;
% opts.colorbar=0;
% opts.title='dark';
% plot_tpm_circle(TPMdark, statenames, opts)
% opts.title='light';
% plot_tpm_circle(TPMlight, statenames, opts)
% 
% % plot_tpm_diff_circle(TPMdark-TPMlight, statenames, 'MinAbsDiff', 0.001);
% % this still looks very messy

    

%% bar graph of dark-light state counts
figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    diffs=results(a).statecounts - results(b).statecounts;
    bh = bar((statenames), diffs);
    bh.FaceColor = 'flat';
    bh.CData = (diffs' <= 0) * [0 0 1] + (diffs' > 0) * [1 0 0];
    ylabel('state count diffs')
    title(sprintf('state count diffs, %s - %s', results(a).condition_name, results(b).condition_name))
end

%% bar graph of dark-light state rates/sec
% to control for different durations in each condition
for c=1:4
    switch results(c).condition_name
        case     'dark laser on'
            results(c).totalframes=sum(dark & laseron);
        case     'light laser on'
            results(c).totalframes=sum(light & laseron);
        case     'dark laser off'
            results(c).totalframes=sum(dark & ~laseron);
        case     'light laser off'
            results(c).totalframes=sum(light & ~laseron);
    end
end

figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
for c=1:4
    a=A(c);
    b=B(c);
    nexttile
    diffs=results(a).statecounts*fps/(results(a).totalframes) - results(b).statecounts*fps/(results(a).totalframes);
    bh = bar((statenames), diffs);
    bh.FaceColor = 'flat';
    bh.CData = (diffs' <= 0) * [0 0 1] + (diffs' > 0) * [1 0 0];
    ylabel('states/second diffs')
    title(sprintf('states/second diffs, %s - %s', results(a).condition_name, results(b).condition_name))
end

%plot bar graph of state rates for each condition
for c=1:4
    figure;
    statespersec=results(c).statecounts*fps/(results(c).totalframes);
    bh = bar((statenames), statespersec);
    bh.FaceColor = 'flat';
    ylabel('states/second')
    title(sprintf('states/second, %s', results(c).condition_name))
end

figure
bh = bar((statenames), [ results(1).statecounts; results(2).statecounts; results(3).statecounts; results(4).statecounts; ]);
ylabel(' state counts ')
title('state counts dark/light/laser')
legend( results(1).condition_name, results(2).condition_name, results(3).condition_name, results(4).condition_name, 'location', 'best');

% %% plot eigenspectrum
% plot_tpm_spectrum(TPMdark)
% 
% %% plot State probability trajectories from each starting state
% 
% %plot for a specific initial condition, such as pause
% pi0=contains(statenames, 'pause'); %0s for all states except 1 at the match
% numsteps=200;
% plot_occupancy_over_time(TPMdark, numsteps,    'stateNames', statenames, 'pi0', pi0);
% 
% %plot for each  initial condition
% for i=1:length(statenames)
%     pi0=zeros(size(statenames));
%     pi0(i)=1;
%     plot_occupancy_over_time(TPMdark, numsteps,    'stateNames', statenames, 'pi0', pi0);
% end
% 











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
