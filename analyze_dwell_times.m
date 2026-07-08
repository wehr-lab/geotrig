%% =========================================================================
%  ANALYZE_DWELL_TIMES
%  Dwell-time (bout duration) analysis per behavioral state, across the
%  four dark/light x laser-on/off condition combinations, plus the four
%  requested difference comparisons:
%    light - dark (laser off)     light - dark (laser on)
%    laser_on - laser_off (dark)  laser_on - laser_off (light)
%
%  Motivation: occupancy (fraction of TIME in a state) can shift even when
%  local drift near that state shows little effect, if the underlying
%  change is in how LONG each visit lasts (dwell time) or how the chain is
%  routed into/out of the state (see the difference-TPM analyses) rather
%  than in local velocity. This was flagged as the natural way to resolve
%  the open question of why dark reduces pause occupancy without a strong
%  local drift signature there.
% =========================================================================

%% 0. ASSUMED INPUTS (frame-level, all Nx1, same length) ---------------------
% state_id (1-6, NaN/0 = NOTA), trial_id, animal_id, dark, laser_on, fs

set(groot, 'defaultAxesFontSize', 18);
statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
dt = 1/fs;

%% =========================================================================
%  1. BUILD BOUT TABLE (run-length encode contiguous same-state segments)
%% =========================================================================
fprintf('=== Building bout table ===\n');

N = numel(state_id);
isChange = [true; (diff(state_id) ~= 0) | (diff(trial_id) ~= 0)];
boutStart = find(isChange);
boutEnd = [boutStart(2:end) - 1; N];
boutState = state_id(boutStart);
boutTrial = trial_id(boutStart);
boutLenFrames = boutEnd - boutStart + 1;

% Exclude NOTA/invalid-state bouts
validState = boutState >= 1 & boutState <= 6 & ~isnan(boutState);

% Exclude left- and right-censored bouts (true duration unknown because
% the bout is cut off by the start/end of its trial's recorded frames).
trial_ids_all = unique(trial_id(~isnan(trial_id)));
firstFrameOfTrial = nan(numel(trial_ids_all),1);
lastFrameOfTrial  = nan(numel(trial_ids_all),1);
for t = 1:numel(trial_ids_all)
    idx = find(trial_id == trial_ids_all(t));
    firstFrameOfTrial(t) = idx(1);
    lastFrameOfTrial(t)  = idx(end);
end
trialFirstMap = containers.Map(trial_ids_all, num2cell(firstFrameOfTrial));
trialLastMap  = containers.Map(trial_ids_all, num2cell(lastFrameOfTrial));

nBouts = numel(boutStart);
isCensored = false(nBouts,1);
nbytes = fprintf('  %.1f%%', 0);
for i = 1:nBouts
    if validState(i) && isKey(trialFirstMap, boutTrial(i))
        if boutStart(i) == trialFirstMap(boutTrial(i)) || boutEnd(i) == trialLastMap(boutTrial(i))
            isCensored(i) = true;
        end
    end
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('  %.1f%%', 100*i/nBouts);
end
fprintf('\n');

keep = validState & ~isCensored;
fprintf('  %d total bouts -> %d after excluding NOTA/invalid (%d) and censored (%d)\n', ...
    nBouts, nnz(keep), nnz(~validState), nnz(isCensored & validState));

DwellSec  = boutLenFrames(keep) * dt;
StateID   = boutState(keep);
TrialID   = boutTrial(keep);
AnimalID  = animal_id(boutStart(keep));
Dark      = logical(dark(boutStart(keep)));
LaserOn   = logical(laser_on(boutStart(keep)));

comboLabels = strings(numel(DwellSec),1);
comboLabels(~Dark & ~LaserOn) = "light_laser_off";
comboLabels(~Dark &  LaserOn) = "light_laser_on";
comboLabels( Dark & ~LaserOn) = "dark_laser_off";
comboLabels( Dark &  LaserOn) = "dark_laser_on";

BoutTbl = table(TrialID, AnimalID, StateID, Dark, LaserOn, comboLabels, DwellSec, ...
    'VariableNames', {'TrialID','AnimalID','StateID','Dark','LaserOn','ConditionCombo','DwellSec'});
BoutTbl.ConditionCombo = categorical(BoutTbl.ConditionCombo, ...
    ["dark_laser_off","dark_laser_on","light_laser_off","light_laser_on"]);

fprintf('  Final bout table: %d bouts across %d trials\n', height(BoutTbl), numel(unique(BoutTbl.TrialID)));

%% =========================================================================
%  2. DESCRIPTIVE PLOT -- dwell time distribution per state, per condition
%% =========================================================================
fprintf('\n=== Plotting dwell-time distributions ===\n');

comboOrder = ["dark_laser_off","dark_laser_on","light_laser_off","light_laser_on"];
figure('Name', 'Dwell time by state and condition combination');
for s = 1:6
    subplot(2,3,s);
    sel = (BoutTbl.StateID == s);
    boxplot(BoutTbl.DwellSec(sel), BoutTbl.ConditionCombo(sel), ...
        'GroupOrder', cellstr(comboOrder), 'Symbol', '');
    set(gca, 'YScale', 'log');
    ylabel('dwell time (s)', 'FontSize', 18);
    title(statenames{s}, 'FontSize', 18);
    grid on;
end
sgtitle('Dwell time by state, across condition combinations (log scale, outliers hidden)', 'FontSize', 18);

%% =========================================================================
%  3. DIFFERENCE TESTS -- the four requested comparisons, per state
%% =========================================================================
fprintf('\n=== Fitting dwell-time difference models ===\n');

light = ~BoutTbl.Dark;   % convenience flag, so contrast coefficients read
                          % directly as "light minus dark" with no sign flip

contrasts = struct( ...
    'name',        {'Light - Dark (laser off)', 'Light - Dark (laser on)', ...
                     'LaserOn - LaserOff (dark)', 'LaserOn - LaserOff (light)'}, ...
    'subsetMask',  {~BoutTbl.LaserOn, BoutTbl.LaserOn, BoutTbl.Dark, light}, ...
    'compareVar',  {light, light, BoutTbl.LaserOn, BoutTbl.LaserOn}, ...
    'compareLabels', {{'dark','light'}, {'dark','light'}, {'laser_off','laser_on'}, {'laser_off','laser_on'}} ...
);

nContrasts = numel(contrasts);
pctChange = nan(6, nContrasts); pctChangeLo = nan(6, nContrasts); pctChangeHi = nan(6, nContrasts);
pValue = nan(6, nContrasts);

for cx = 1:nContrasts
    fprintf('\n--- %s ---\n', contrasts(cx).name);
    betaName = sprintf('CompareVar_%s', contrasts(cx).compareLabels{2});
    for s = 1:6
        sel = (BoutTbl.StateID == s) & contrasts(cx).subsetMask;
        subTbl = BoutTbl(sel,:);
        subTbl.CompareVar = categorical(contrasts(cx).compareVar(sel), [false true], contrasts(cx).compareLabels);
        subTbl.LogDwell = log(subTbl.DwellSec);

        if nnz(subTbl.CompareVar == contrasts(cx).compareLabels{1}) < 20 || ...
           nnz(subTbl.CompareVar == contrasts(cx).compareLabels{2}) < 20
            fprintf('  %-12s insufficient bouts, skipped\n', statenames{s});
            continue
        end

        lme = fitlme(subTbl, 'LogDwell ~ CompareVar + (1|AnimalID) + (1|TrialID)');
        rowIdx = strcmp(lme.Coefficients.Name, betaName);
        c = lme.Coefficients(rowIdx,:);

        % Convert log-scale beta (and CI) to a percent change in dwell time
        pctChange(s,cx)   = (exp(c.Estimate) - 1) * 100;
        pctChangeLo(s,cx) = (exp(c.Estimate - 1.96*c.SE) - 1) * 100;
        pctChangeHi(s,cx) = (exp(c.Estimate + 1.96*c.SE) - 1) * 100;
        pValue(s,cx) = c.pValue;

        fprintf('  %-12s %+6.1f%% [%+6.1f%%, %+6.1f%%], p=%.4g (n=%d bouts)\n', ...
            statenames{s}, pctChange(s,cx), pctChangeLo(s,cx), pctChangeHi(s,cx), pValue(s,cx), height(subTbl));
    end
end

%% =========================================================================
%  4. DIFFERENCE SUMMARY PLOT
%% =========================================================================
figure('Name', 'Dwell-time percent change by contrast');
for cx = 1:nContrasts
    subplot(2,2,cx); hold on;
    err = [pctChange(:,cx) - pctChangeLo(:,cx), pctChangeHi(:,cx) - pctChange(:,cx)];
    b = bar(1:6, pctChange(:,cx));
    errorbar(1:6, pctChange(:,cx), err(:,1), err(:,2), 'k.', 'LineWidth', 1.5, 'CapSize', 8);
    yline(0, '--k', 'HandleVisibility','off');
    set(gca, 'XTick', 1:6, 'XTickLabel', statenames, 'XTickLabelRotation', 30);
    ylabel('% change in dwell time', 'FontSize', 18);
    title(contrasts(cx).name, 'FontSize', 18);
    grid on;
end
sgtitle('Dwell-time percent change, with 95% CI, across the four contrasts', 'FontSize', 18);

fprintf('\nDone.\n');

