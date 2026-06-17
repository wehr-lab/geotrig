% Decision Gate Framework

% Assuming total frames N
N=num_geoframes;
% 1. Create a label vector 'state_vec' (1=FastChase, 2=MedChase, ..., 6=Wander)

statenames={ ...
    'hot pursuit' , ...
    'chase' , ...
    'follow' , ...
    'stalk' , ...
    'pause' , ...
    'wander' };


state_vec = NaN(N, 1);
state_vec(hotpursuit) = 1;
state_vec(chase) = 2;
state_vec(follow) = 3;
state_vec(stalk) = 4;
state_vec(pause) = 5;
state_vec(wander) = 6;


% 2. Create an event vector 'event_vec'


eventnames = {'failed_approach','contact_loss','intercept','cricket_jump', 'rangemin'};


event_vec = zeros(N, 1);
event_vec(failed_approach_event_frames) = 1; %
event_vec(contact_loss_event_frames) = 2; %
event_vec(intercept_event_frames) = 3; %
event_vec(cricket_jump_event_frames) = 4; %
event_vec(rangemin_event_frames) = 5; %




%% statistics: test whether event affects persistence N steps ahead, using glm

if 0 
for event_type=1:eventnames

    % Identify frames where 'Lose Contact' happened during a chase
    %    chase_states = [1, 2, 3]; % Fast, Med, Slow
    chase_states = [4:6]; % Fast, Med, Slow

    event_indices = find(event_vec == event_type);

    fprintf('\n____________________________________________')
    fprintf('\nmeasuring chase states persistence following event %s', eventnames{event_type})

    persistence = []; % 0 for reset, 1 for persist
    cond_list = {};
    for i = 1:length(event_indices)
        idx = event_indices(i);
        if ismember(state_vec(idx), chase_states)
            % Look at state 10 frames after the event
            if ismember(state_vec(idx + 200*2), chase_states)
                persistence(end+1) = 1; % Persisted
            else
                persistence(end+1) = 0; % Reset
            end
            % Grab the condition at that frame
            if dark(idx) & laseron(idx)
                cond_list{end+1} = 'dark-laserON';
            elseif dark(idx) & ~laseron(idx)
                cond_list{end+1} = 'dark-laserOFF';
            elseif light(idx) & laseron(idx)
                cond_list{end+1} = 'light-laserON';
            elseif light(idx) & ~laseron(idx)
                cond_list{end+1} = 'light-laserOFF';
            end

        end
    end


    % 2. Now you have a list matching your event_indices length
    Condition = categorical(cond_list);

    % 3. Build your final table for the GLM
    T = table(Condition', persistence', 'VariableNames', {'Condition', 'Persistence'});

    % Set 'light-laserOFF' as the reference
    T.Condition = reordercats(T.Condition, { 'light-laserOFF', 'dark-laserOFF', 'dark-laserON', 'light-laserON'});

    fprintf('\nmean persistence in each experimental condition group:')
    groupsummary(T, 'Condition', 'mean', 'Persistence')

    % You can now run a Logistic Regression:
    model = fitglm(T, 'Persistence ~ Condition', 'Distribution', 'binomial')

end

% 1. Interpreting the Estimates (Log-Odds)The Intercept (dark-laserOFF):
% $4.8752$. This is a very high positive value, meaning that in the
% dark-laserOFF condition, the probability of the mouse persisting is
% extremely high (the log-odds are strongly positive).Condition Estimates:
% All other estimates (0.146, 0.879, 0.144) are relatively small when
% compared to the standard error (SE).The "Chi-Square" Result: The most
% important line in your output is the bottom one: "Chi^2-statistic vs.
% constant model: 1.27, p-value = 0.737." * This means your experimental
% conditions (Light/Dark and Laser) are NOT significantly affecting
% persistence. Your model is no better at predicting persistence than a
% model that simply predicts the average persistence for every trial
% regardless of condition.2. What this means for your BiologyCeiling
% Effect: Since the intercept is so high ($4.8752$), your mice are likely
% "maxed out" on persistence in the baseline condition. If they are already
% persisting $95\%+$ of the time, it is statistically impossible for any
% manipulation to significantly increase that behavior further.Lack of
% Effect: Based on the high p-values ($0.84, 0.31, 0.85$), you cannot
% reject the null hypothesis. Lighting and Laser, in your current
% experimental setup, do not appear to be driving the decision to "reset"
% vs "stay."

% Check for a "Ceiling Effect": Calculate the mean persistence probability
% for each group. If they are all hovering around $0.95$ to $0.99$, your
% experiment is likely hitting a behavioral limit where the mouse always
% persists, regardless of external conditions.

% Inspect the "Reset" Events: If the mouse almost never resets, your
% persistence vector is almost all 1s. Statistical models struggle when
% there is no variance (i.e., when the outcome is essentially a constant).

% Narrow the Window: Are you defining persistence too broadly? If "persist"
% means "stayed within 10 seconds," perhaps the mouse is always staying
% within 10 seconds. Try shortening the time window (e.g., look for
% persistence within 1 second) to see if you can find a condition where the
% mouse does reset.

% positive beta (Estimate) means the condition increases the log-odds of persistence
% negative beta (Estimate) means the condition decreases the log-odds of persistence
% 0 beta means no effect of condition on persistence

% Predict probabilities for each group
T_levels = table(categorical({'dark-laserOFF', 'dark-laserON', 'light-laserOFF', 'light-laserON'}'), ...
    'VariableNames', {'Condition'});
[y_pred, y_ci] = predict(model, T_levels);

% Plot as a grouped interaction
figure
bar(y_pred);
hold on;
errorbar(1:4, y_pred, y_pred-y_ci(:,1), y_ci(:,2)-y_pred, 'k', 'LineStyle', 'none');
set(gca, 'XTickLabel', {'Dark-Off', 'Dark-On', 'Light-Off', 'Light-On'});
ylabel('Probability of Persistence');
title('Effect of Lighting and Laser on Mouse Persistence');

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot effect of event on persistence


mean_persistence = nan(length(eventnames), length(window)); % 0 for reset, 1 for persist
mean_persistence_dark_laserON=mean_persistence;
mean_persistence_dark_laserOFF=mean_persistence;
mean_persistence_light_laserON=mean_persistence;
mean_persistence_light_laserOFF=mean_persistence;
win_start=-1; %in seconds
win_stop=5;%in seconds
window=win_start*200:win_stop*200;
chase_states = [1, 2, 3]; % Fast, Med, Slow
search_states = [4, 5, 6]; % search

init_states=chase_states; %condition on event occurring in one of these states
init_states_label='chase';
persist_states=chase_states; %define persistence as being in one of these states
persist_states_label='chase';
% init_states=search_states; %condition on event occurring in one of these states
% init_states_label='search';
% persist_states=search_states; %define persistence as being in one of these states
% persist_states_label='search';

for event_type=1:length(eventnames)

    % Identify frames where 'Lose Contact' happened during a chase
    event_indices = find(event_vec == event_type);

    fprintf('\n____________________________________________')
    fprintf('\nmeasuring chase states persistence following event %s', eventnames{event_type})

    persistence=nan(length(event_indices), length(window));
    persistence_dark_laserON=nan(size(window));
    persistence_dark_laserOFF=nan(size(window));
    persistence_light_laserON=nan(size(window));
    persistence_light_laserOFF=nan(size(window));
    cond_list = {};
    for i = 1:length(event_indices)
        idx = event_indices(i);
        if ismember(state_vec(idx), init_states)

            persistence(i, 1:length(window))= ismember(state_vec(idx+window), persist_states);

            % Grab the condition at that frame
            if dark(idx) & laseron(idx)
                persistence_dark_laserON(size(persistence_dark_laserON, 1)+1, 1:length(window))= ismember(state_vec(idx+window), persist_states);
            elseif dark(idx) & ~laseron(idx)
                persistence_dark_laserOFF(size(persistence_dark_laserOFF, 1)+1, 1:length(window))= ismember(state_vec(idx+window), persist_states);
            elseif light(idx) & laseron(idx)
                persistence_light_laserON(size(persistence_light_laserON, 1)+1, 1:length(window))= ismember(state_vec(idx+window), persist_states);
            elseif light(idx) & ~laseron(idx)
                persistence_light_laserOFF(size(persistence_light_laserOFF, 1)+1, 1:length(window))= ismember(state_vec(idx+window), persist_states);
            end

        end
    end

    mean_persistence(event_type, :)=smooth(nanmean(persistence), 100);
    mean_persistence_dark_laserON(event_type, :)=smooth(nanmean(persistence_dark_laserON), 100);
    mean_persistence_dark_laserOFF(event_type, :)=smooth(nanmean(persistence_dark_laserOFF), 100);
    mean_persistence_light_laserON(event_type, :)=smooth(nanmean(persistence_light_laserON), 100);
    mean_persistence_light_laserOFF(event_type, :)=smooth(nanmean(persistence_light_laserOFF), 100);
end



figure
plot(window/200, mean_persistence,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, all conditions', persist_states_label))


figure
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 2*420])

nexttile
plot(window/200, mean_persistence_dark_laserON,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, dark_laserON', persist_states_label), 'interpreter', 'none')

nexttile
plot(window/200, mean_persistence_dark_laserOFF,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, dark_laserOFF', persist_states_label), 'interpreter', 'none')

nexttile
plot(window/200, mean_persistence_light_laserON,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, light_laserON', persist_states_label), 'interpreter', 'none')

nexttile
plot(window/200, mean_persistence_light_laserOFF,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, light_laserOFF', persist_states_label), 'interpreter', 'none')


%%%%%%%%%%
%differences

figure
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 2*420])
yl=[-.4 .4];


nexttile
plot(window/200, mean_persistence_dark_laserON-mean_persistence_dark_laserOFF,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['difference persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, dark_laserON-dark_laserOFF', persist_states_label), 'interpreter', 'none')
ylim(yl)

nexttile
plot(window/200, mean_persistence_light_laserON-mean_persistence_light_laserOFF,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['difference persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, light_laserON-light_laserOFF', persist_states_label), 'interpreter', 'none')
ylim(yl)

nexttile
plot(window/200, mean_persistence_dark_laserOFF-mean_persistence_light_laserOFF,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['difference persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, dark_laserOFF-light_laserOFF', persist_states_label), 'interpreter', 'none')
ylim(yl)

nexttile
plot(window/200, mean_persistence_dark_laserON-mean_persistence_light_laserON,  'linew', 2)
xticks(win_start:win_stop)
grid on
xlabel(['time relative to event in ', init_states_label])
ylabel(['difference persistence in ', persist_states_label])
legend(eventnames, 'interpreter', 'none')
set(gca, 'fontsize', 18)
title(sprintf('Persistence in %s, dark_laserON-light_laserON', persist_states_label), 'interpreter', 'none')
ylim(yl)

