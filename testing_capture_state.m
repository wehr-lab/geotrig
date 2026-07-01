%what does it look like if we add capture as a state ?

capture=false(size(chase));
capture(capture_frames)=true;

stateMask = [    hotpursuit(:),chase(:),follow(:), stalk(:), wander(:), pause(:), capture(:), noneoftheabove(:)];

stateMask(117153+[-10:10],:)

win=[115410-10:117153+10];

figure
hold on
plot(win, az(win), 'm')
plot(win, mouse_spd(win), 'r')
plot(win, range(win), 'g')


noneIdx=find(contains(statenames,'none'));
captureIdx=find(contains(statenames,'capture'));
M = zeros(length(find(state_sequence==captureIdx)), 11+2); %
j=0;

for i = 1:length(state_sequence) - 1
    current_state = state_sequence(i);
    next_state = state_sequence(i+1);
    if current_state == noneIdx || next_state == noneIdx %skip counting none transitions
        %continue;
    end
    if current_state == captureIdx  %skip counting transitions from capture to anything
        %(because capture terminates a trial)
        % fprintf('.')
        j=j+1;
        M(j,:)=state_sequence(i-10:i+2);
    end

end

%%
 statenames={ ...
        'hot pursuit' , ...
        'chase' , ...
        'following' , ...
        'stalk' , ...
        'wander' , ...
        'pause' , ...
        'capture', ...
        'none'};
    num_states = length(statenames); % Chase, Pause, Wander, etc.

T = zeros(num_states, num_states); % The raw count matrix

noneIdx=find(contains(statenames,'none'));
captureIdx=find(contains(statenames,'capture'));

% Count transitions
for i = 1:length(state_sequence) - 1
    current_state = state_sequence(i);
    next_state = state_sequence(i+1);
    % if current_state == noneIdx || next_state == noneIdx %skip counting none transitions
    %     %continue;
    % end
    if current_state == captureIdx  %skip counting transitions from capture to anything
        continue
    end
    if next_state == captureIdx
        if current_state==noneIdx
            current_state = state_sequence(i-1);
            %apparently almost every capture frame is preceded by "none of the above",
            %perhaps because things get messy right at the end
            %so here we skip over "none" to find the previous state
        end
    end

    T(current_state, next_state) = T(current_state, next_state) + 1;
end

keepidx=setdiff(1:num_states, noneIdx);
T=T(keepidx, keepidx);
statenames=statenames(keepidx);
states=states(keepidx);
num_states = length(states);


figure;
hT=heatmap(statenames, statenames, T, ...
    'Title', ['Behavioral Transition Counts, ', condition_name], ...
    'Colormap', parula, 'GridVisible', 'off');
xlabel('To State')
ylabel('From State')

T = T + 0.01;

    TPM = T ./ sum(T, 2);

TPM(captureIdx,:)=0;

    figure
    

    hTPM=heatmap(statenames, statenames, TPM, ...
        'Title', ['Behavioral Transition Probabilities (rows sum to 1), ', condition_name], ...
        'Colormap', parula, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
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

%plot without clustering
figure
hTPM_uc=heatmap(statenames, statenames, TPM, ...
    'Title', ['TPM dark, no clustering'], ...
    'Colormap', parula, 'GridVisible', 'off', 'CellLabelFormat', '%0.2g', 'CellLabelColor', 'none'); %'auto'
xlabel('To State')
ylabel('From State')
colormap hot

%% plot tpm circle
opts.minProb=0.1;
opts.colorbar=0;
opts.title='dark';
plot_tpm_circle(TPM, statenames, opts)
