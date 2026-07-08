%% =========================================================================
%  FIT_CONDITION_AND_EVENT_EFFECTS
%  OU-process modeling of z(t) for condition comparisons (light/dark,
%  laser on/off), plus peri-event trigger analysis and jump regression
%  for point-process events (cricket jump, intercept, etc.)
%
%  Model: dz = -theta*(z - mu)*dt + sigma*dW      (Ornstein-Uhlenbeck)
%    theta = relaxation rate (how fast z returns to baseline)
%    mu    = equilibrium engagement level
%    sigma = volatility along the axis
%
%  Fitting unit is the TRIAL/BOUT, not the frame -- this respects
%  within-trial autocorrelation and lets condition/animal be proper
%  grouping factors in a mixed model, rather than pooling autocorrelated
%  frames as if independent.
% =========================================================================

%% 0. ASSUMED INPUTS ---------------------------------------------------------
% z          : Nx1, latent engagement trajectory from estimate_engagement_axis.m
%              (NaN where invalid/excluded)
% trial_id   : Nx1, integer/categorical ID grouping frames into contiguous
%              bouts. z(t) should NOT be treated as continuous across a
%              trial boundary (different recording, animal, or long gap).
% condition  : Nx1 (or table indexed by trial), e.g. categorical
%              {'light','dark'} or {'laser_on','laser_off'}
% animal_id  : Nx1 (or per-trial), grouping factor for mixed models
% fs         : scalar, frame rate (Hz) -- you said 200 fps elsewhere
%
% event_tbl  : table with one row per point-process event, columns:
%              EventFrame (index into z/trial_id), EventType (categorical:
%              'failed_approach','contact_loss','contact_gain','intercept',
%              'cricket_jump','range_min'), TrialID, AnimalID

fprintf('=== [0/5] Setup ===\n');
fs=fps;
dt = 1/fs;
fprintf('  dt = %.5f s (fs = %.1f Hz)\n', dt, fs);

if ~exist('animal_id')
    derive_trial_and_animal_id
end
if ~exist('event_tbl')
    build_event_table
end
laser_on=laseron;
laser_off=~laseron;

condition=categorical(zeros(size(az)), [0], {'none'});
condition(dark & laser_off)={'dark & laser_off'};
condition(dark & laser_on)={'dark & laser_on'};
condition(light & laser_off)={'light & laser_off'};
condition(light & laser_on)={'light & laser_on'};

fprintf('%d frames still uncategorized (should be 0)\n', sum(condition == 'none'));
if sum(condition == 'none')>0
    error (sprintf('%d frames still uncategorized (should be 0)\n', sum(condition == 'none')))
end
condition = removecats(condition, 'none');   % drop it once confirmed empty

%% =========================================================================
%  PART A -- OU FITTING PER TRIAL, CONDITION COMPARISON
%% =========================================================================

isDark = TrialParams.Condition == 'dark & laser_off' | TrialParams.Condition == 'dark & laser_on';
darkOnly = TrialParams(isDark, :);
darkOnly.Condition = removecats(darkOnly.Condition);   % drop the 2 unused light levels

lme_laser_in_dark = fitlme(darkOnly, 'Theta ~ Condition + (1|AnimalID)');
disp(lme_laser_in_dark.Coefficients);

% in the A1 loop, alongside the existing Cond/Animal extraction:
darkVals = dark(idx); darkVals = darkVals(~isnan(double(darkVals)));  % or just dark(idx)
if ~isempty(darkVals), Dark(i) = darkVals(1); end
laserVals = laser_on(idx);
if ~isempty(laserVals), LaserOn(i) = laserVals(1); end

darkOnly = TrialParams(TrialParams.Dark, :);
lme_laser_in_dark = fitlme(darkOnly, 'Theta ~ LaserOn + (1|AnimalID)');

%% A1. Fit OU parameters within each trial via discretized AR(1) --------
fprintf('\n=== [1/5] Part A1: fitting OU parameters per trial ===\n');
trial_ids = unique(trial_id(~isnan(trial_id)));
minFramesPerTrial = 30;   % need enough frames for a stable phi/theta estimate

nT = numel(trial_ids);
Theta = nan(nT,1); Mu = nan(nT,1); Sigma2 = nan(nT,1);
R2 = nan(nT,1); Nframes = nan(nT,1);
Cond = strings(nT,1); Animal = nan(nT,1);

nbytes = fprintf('  %.1f%%', 0);
for i = 1:nT
    idx = (trial_id == trial_ids(i));
    zseg = z(idx);
    if nnz(~isnan(zseg)) >= minFramesPerTrial
        fit = local_fitOU_AR1(zseg, dt);
        Theta(i) = fit.theta; Mu(i) = fit.mu; Sigma2(i) = fit.sigma2;
        R2(i) = fit.R2; Nframes(i) = fit.n;

        condVals = condition(idx); condVals = condVals(~ismissing(condVals));
        if ~isempty(condVals), Cond(i) = string(condVals(1)); end
        animVals = animal_id(idx); animVals = animVals(~isnan(animVals));
        if ~isempty(animVals), Animal(i) = animVals(1); end
    end
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('  %.1f%%', 100*i/nT);
end
fprintf('\n');

TrialParams = table(trial_ids, Animal, Cond, Theta, Mu, Sigma2, R2, Nframes, ...
    'VariableNames', {'TrialID','AnimalID','Condition','Theta','Mu','Sigma2','R2','N'});

% Drop trials with unstable/non-mean-reverting fits (phi outside (0,1))
bad = isnan(TrialParams.Theta) | TrialParams.Theta <= 0 | ~isfinite(TrialParams.Theta);
fprintf('  Dropped %d / %d trials (too short or non-mean-reverting fit).\n', nnz(bad), nT);
TrialParams(bad,:) = [];
TrialParams.Condition = categorical(TrialParams.Condition);

%% A2. Condition comparison -----------------------------------------------
fprintf('\n=== [2/5] Part A2: condition comparison ===\n');
haveStats = license('test','statistics_toolbox') && exist('fitlme','file') == 2;

if haveStats
    fprintf('  Fitting mixed models (Statistics Toolbox) ...\n');
    lme_theta = fitlme(TrialParams, 'Theta ~ Condition + (1|AnimalID)');
    lme_mu    = fitlme(TrialParams, 'Mu ~ Condition + (1|AnimalID)');
    TrialParams.LogSigma2 = log(TrialParams.Sigma2);
    lme_sig   = fitlme(TrialParams, 'LogSigma2 ~ Condition + (1|AnimalID)');

    disp('theta ~ condition:'); disp(lme_theta.Coefficients);
    disp('mu ~ condition:');    disp(lme_mu.Coefficients);
    disp('log(sigma^2) ~ condition:'); disp(lme_sig.Coefficients);
else
    fprintf('  Statistics Toolbox / fitlme not found -- using cluster bootstrap fallback.\n');
    nboot = 2000;
    for pname = ["Theta","Mu","Sigma2"]
        fprintf('  Bootstrapping %s (%d resamples) ...\n', pname, nboot);
        [diffEst, ci] = local_clusterBootstrapDiff(TrialParams, pname, nboot);
        fprintf('  %-8s condition diff: %+.4f  95%% CI [%+.4f, %+.4f]\n', ...
            pname, diffEst, ci(1), ci(2));
    end
end

%% =========================================================================
%  PART B -- EVENT-TRIGGERED ANALYSIS (jump-diffusion / point-process link)
%% =========================================================================

%% B1. Peri-event traces (PSTH-style, for continuous z) --------------------
fprintf('\n=== [3/5] Part B1: peri-event traces ===\n');
nPre  = round(1.0*fs);   % 1 s before event
nPost = round(1.0*fs);   % 1 s after event
tAxis = (-nPre:nPost) / fs;

eventTypes = categories(categorical(event_tbl.EventType));
figure; hold on;
traceByType = struct();
for k = 1:numel(eventTypes)
    et = eventTypes{k};
    fprintf('  Aligning event type %d/%d: %s ...\n', k, numel(eventTypes), et);
    rows = strcmp(string(event_tbl.EventType), et);
    frames = event_tbl.EventFrame(rows);
    evTrial = event_tbl.TrialID(rows);
    M = local_periEventTrace(z, trial_id, frames, evTrial, nPre, nPost);
    traceByType.(matlab.lang.makeValidName(et)) = M;
    mtrace = mean(M, 1, 'omitnan');
    setrace = std(M, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(M),1));
    plot(tAxis, mtrace, 'DisplayName', sprintf('%s (n=%d)', et, size(M,1)), 'linewidth', 2);
    fill([tAxis fliplr(tAxis)], [mtrace+setrace fliplr(mtrace-setrace)], 'k', ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility','off');
end
xline(0, '--k', 'HandleVisibility','off');
xlabel('Time from event (s)'); ylabel('z(t), engagement axis');
title('Peri-event latent trajectory by event type');
legend('Location','bestoutside', 'interpreter', 'none'); grid on;
set(gca, 'fontsize', 18)
set(gcf, 'pos', [     451         648        1109         590]);

%% B2. Event-triggered jump size (Delta z), pre vs post window -------------
fprintf('\n=== [4/5] Part B2: event-triggered jump size (Delta z) ===\n');
preWin  = [-round(0.5*fs), -round(0.15*fs)];  % frames relative to event
postWin = [ round(0.15*fs),  round(0.5*fs)];  % gap around 0 avoids the
                                                % event itself contaminating
                                                % both windows
fprintf('  Computing Delta z for %d events ...\n', height(event_tbl));
EventEffects = local_deltaZfromEvents(z, trial_id, event_tbl, preWin, postWin);

%% B3. Test each event type's Delta z against zero --------------------------
fprintf('\n=== [5/5] Part B3: testing event-triggered jump size ===\n');
if haveStats
    fprintf('  Fitting mixed model (Statistics Toolbox) ...\n');
    EventEffects.EventType = categorical(EventEffects.EventType);
    lme_event = fitlme(EventEffects, 'DeltaZ ~ EventType + (1|AnimalID)');
    disp(lme_event.Coefficients);
    disp(anova(lme_event));
else
    fprintf('  Statistics Toolbox / fitlme not found -- using permutation test fallback.\n');
    for k = 1:numel(eventTypes)
        et = eventTypes{k};
        rows = strcmp(string(EventEffects.EventType), et);
        dz = EventEffects.DeltaZ(rows);
        dz = dz(~isnan(dz));
        if numel(dz) < 5, continue; end
        fprintf('  Permutation testing event type %d/%d: %s ...\n', k, numel(eventTypes), et);
        [pNull, nullMean] = local_permTestEventEffect(z, trial_id, event_tbl, et, ...
            preWin, postWin, 1000);
        fprintf('  %-16s mean dz = %+.4f (n=%d), null mean = %+.4f, p = %.4f\n', ...
            et, mean(dz), numel(dz), nullMean, pNull);
    end
end

fprintf('\nDone.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function out = local_fitOU_AR1(zseg, dt)
    % Discretized OU (AR1) fit: z(t+1) = phi*z(t) + c + eps
    X = zseg(1:end-1); Y = zseg(2:end);
    valid = ~isnan(X) & ~isnan(Y);
    n = nnz(valid);
    A = [X(valid), ones(n,1)];
    coefs = A \ Y(valid);
    phi = coefs(1); c = coefs(2);
    resid = Y(valid) - A*coefs;
    residVar = var(resid);

    if phi <= 0 || phi >= 1 || ~isfinite(phi)
        out = struct('theta',NaN,'mu',NaN,'sigma2',NaN,'R2',NaN,'n',n);
        return
    end
    theta = -log(phi) / dt;
    mu = c / (1 - phi);
    sigma2 = residVar * 2*theta / (1 - phi^2);
    R2 = 1 - residVar / var(Y(valid));
    out = struct('theta',theta,'mu',mu,'sigma2',sigma2,'R2',R2,'n',n);
end

function [diffEst, ci] = local_clusterBootstrapDiff(T, pname, nboot)
    % Two-stage cluster bootstrap: resample animals with replacement,
    % average trial-level parameter within each resampled animal, then
    % compute the between-condition mean difference. Repeat nboot times.
    conds = categories(T.Condition);
    assert(numel(conds) == 2, 'This fallback assumes exactly 2 conditions.');
    animals = unique(T.AnimalID);
    diffs = nan(nboot,1);
    nbytes = fprintf('    %.1f%%', 0);
    for b = 1:nboot
        samp = animals(randi(numel(animals), numel(animals), 1));
        m = nan(numel(samp), 2);
        for a = 1:numel(samp)
            for cix = 1:2
                rows = T.AnimalID == samp(a) & T.Condition == conds{cix};
                m(a,cix) = mean(T.(pname)(rows), 'omitnan');
            end
        end
        diffs(b) = mean(m(:,1) - m(:,2), 'omitnan');
        fprintf(repmat('\b',1,nbytes));
        nbytes = fprintf('    %.1f%%', 100*b/nboot);
    end
    fprintf('\n');
    diffEst = mean(diffs, 'omitnan');
    ci = prctile(diffs, [2.5, 97.5]);
end

function M = local_periEventTrace(z, trial_id, eventFrames, eventTrialID, nPre, nPost)
    % Aligns z around each event, padding with NaN if the window would
    % cross a trial boundary or go out of bounds.
    n = numel(eventFrames);
    L = nPre + nPost + 1;
    M = nan(n, L);
    showProgress = n > 200;
    if showProgress, nbytes = fprintf('    %.1f%%', 0); end
    for i = 1:n
        f = eventFrames(i);
        idx = (f-nPre):(f+nPost);
        inBounds = idx >= 1 & idx <= numel(z);
        row = nan(1, L);
        validIdx = idx(inBounds);
        sameTrial = trial_id(validIdx) == eventTrialID(i);
        rowVals = nan(1, numel(validIdx));
        rowVals(sameTrial) = z(validIdx(sameTrial));
        row(inBounds) = rowVals;
        M(i,:) = row;
        if showProgress
            fprintf(repmat('\b',1,nbytes));
            nbytes = fprintf('    %.1f%%', 100*i/n);
        end
    end
    if showProgress, fprintf('\n'); end
end

function EventEffects = local_deltaZfromEvents(z, trial_id, event_tbl, preWin, postWin)
    n = height(event_tbl);
    DeltaZ = nan(n,1); PreMean = nan(n,1); PostMean = nan(n,1);
    showProgress = n > 200;
    if showProgress, nbytes = fprintf('    %.1f%%', 0); end
    for i = 1:n
        f = event_tbl.EventFrame(i);
        tid = event_tbl.TrialID(i);
        preIdx  = (f+preWin(1)):(f+preWin(2));
        postIdx = (f+postWin(1)):(f+postWin(2));
        if ~(any(preIdx < 1) || any(postIdx > numel(z))) && ...
           ~(any(trial_id(preIdx) ~= tid) || any(trial_id(postIdx) ~= tid))
            PreMean(i)  = mean(z(preIdx), 'omitnan');
            PostMean(i) = mean(z(postIdx), 'omitnan');
            DeltaZ(i)   = PostMean(i) - PreMean(i);
        end
        if showProgress
            fprintf(repmat('\b',1,nbytes));
            nbytes = fprintf('    %.1f%%', 100*i/n);
        end
    end
    if showProgress, fprintf('\n'); end
    EventEffects = table(event_tbl.TrialID, event_tbl.AnimalID, ...
        string(event_tbl.EventType), DeltaZ, PreMean, PostMean, ...
        'VariableNames', {'TrialID','AnimalID','EventType','DeltaZ','PreMean','PostMean'});
end

function [pval, nullMean] = local_permTestEventEffect(z, trial_id, event_tbl, eventType, preWin, postWin, nperm)
    % Null: recompute Delta z at random frames within the same trials that
    % actually contain this event type, instead of at the real event times.
    rows = strcmp(string(event_tbl.EventType), eventType);
    realTbl = event_tbl(rows,:);
    obsDz = local_deltaZfromEvents(z, trial_id, realTbl, preWin, postWin);
    obsMean = mean(obsDz.DeltaZ, 'omitnan');

    margin = max(abs([preWin postWin])) + 1;
    nullMeans = nan(nperm,1);
    nbytes = fprintf('    %.1f%%', 0);
    for p = 1:nperm
        permTbl = realTbl;
        for i = 1:height(permTbl)
            tid = permTbl.TrialID(i);
            candidates = find(trial_id == tid);
            candidates = candidates(candidates > margin & candidates < numel(z)-margin);
            if isempty(candidates), continue; end
            permTbl.EventFrame(i) = candidates(randi(numel(candidates)));
        end
        permDz = local_deltaZfromEvents(z, trial_id, permTbl, preWin, postWin);
        nullMeans(p) = mean(permDz.DeltaZ, 'omitnan');
        fprintf(repmat('\b',1,nbytes));
        nbytes = fprintf('    %.1f%%', 100*p/nperm);
    end
    fprintf('\n');
    nullMean = mean(nullMeans, 'omitnan');
    pval = mean(abs(nullMeans - nullMean) >= abs(obsMean - nullMean));
end