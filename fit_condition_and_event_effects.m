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
% dark       : Nx1 logical, true when the frame is in the dark condition
% laser_on   : Nx1 logical, true when the frame is in the laser-on condition
%              (dark & laser_on are independent binary factors -- this
%              script fits their interaction rather than collapsing them
%              into a single 4-level combined category)
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


%% =========================================================================
%  PART A -- OU FITTING PER TRIAL, CONDITION COMPARISON
%% =========================================================================

%% A1. Fit OU parameters within each trial via discretized AR(1) --------
fprintf('\n=== [1/5] Part A1: fitting OU parameters per trial ===\n');
trial_ids = unique(trial_id(~isnan(trial_id)));
minFramesPerTrial = 30;   % need enough frames for a stable phi/theta estimate

nT = numel(trial_ids);
Theta = nan(nT,1); Mu = nan(nT,1); Sigma2 = nan(nT,1);
R2 = nan(nT,1); Nframes = nan(nT,1);
Dark = false(nT,1); LaserOn = false(nT,1); Animal = nan(nT,1);

nbytes = fprintf('  %.1f%%', 0);
for i = 1:nT
    idx = (trial_id == trial_ids(i));
    zseg = z(idx);
    if nnz(~isnan(zseg)) >= minFramesPerTrial
        fit = local_fitOU_AR1(zseg, dt);
        Theta(i) = fit.theta; Mu(i) = fit.mu; Sigma2(i) = fit.sigma2;
        R2(i) = fit.R2; Nframes(i) = fit.n;

        darkVals = dark(idx); darkVals = darkVals(~isnan(double(darkVals)));
        if ~isempty(darkVals), Dark(i) = darkVals(1); end
        laserVals = laser_on(idx); laserVals = laserVals(~isnan(double(laserVals)));
        if ~isempty(laserVals), LaserOn(i) = laserVals(1); end
        if numel(unique(darkVals)) > 1 || numel(unique(laserVals)) > 1
            warning('Trial %d has a condition that changes mid-trial -- using first frame''s value.', trial_ids(i));
        end
        animVals = animal_id(idx); animVals = animVals(~isnan(animVals));
        if ~isempty(animVals), Animal(i) = animVals(1); end
    end
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('  %.1f%%', 100*i/nT);
end
fprintf('\n');

TrialParams = table(trial_ids, Animal, Dark, LaserOn, Theta, Mu, Sigma2, R2, Nframes, ...
    'VariableNames', {'TrialID','AnimalID','Dark','LaserOn','Theta','Mu','Sigma2','R2','N'});

% Drop trials with unstable/non-mean-reverting fits (phi outside (0,1))
bad = isnan(TrialParams.Theta) | TrialParams.Theta <= 0 | ~isfinite(TrialParams.Theta);
fprintf('  Dropped %d / %d trials (too short or non-mean-reverting fit).\n', nnz(bad), nT);
TrialParams(bad,:) = [];

% Set reference levels explicitly (dark=false / laser_off is the reference
% for both factors) so coefficients below have a fixed, known meaning.
TrialParams.Dark    = categorical(TrialParams.Dark,    [false true], {'light','dark'});
TrialParams.LaserOn = categorical(TrialParams.LaserOn, [false true], {'laser_off','laser_on'});

%% A2. Condition comparison -----------------------------------------------
fprintf('\n=== [2/5] Part A2: condition comparison ===\n');

% --- Full 2x2 interaction model (uses all trials at once) ------------------
% With Dark='light' and LaserOn='laser_off' as reference levels:
%   Intercept                   -> light, laser_off baseline
%   Dark_dark                    -> dark vs light, WHEN laser is off
%   LaserOn_laser_on              -> laser effect, WHEN in the light
%   Dark_dark:LaserOn_laser_on     -> how much the laser effect CHANGES in the dark
% i.e. the laser effect *in the dark* is NOT read directly off this table --
% it's (LaserOn_laser_on + Dark_dark:LaserOn_laser_on). See the subset
% comparisons below for that number directly, with no arithmetic required.
fprintf('  Fitting Theta/Mu/Sigma2 ~ Dark*LaserOn + (1|AnimalID) ...\n');
lme_theta = fitlme(TrialParams, 'Theta ~ Dark*LaserOn + (1|AnimalID)');
lme_mu    = fitlme(TrialParams, 'Mu ~ Dark*LaserOn + (1|AnimalID)');
TrialParams.LogSigma2 = log(TrialParams.Sigma2);
lme_sig   = fitlme(TrialParams, 'LogSigma2 ~ Dark*LaserOn + (1|AnimalID)');

disp('theta ~ dark*laser:'); disp(lme_theta.Coefficients);
disp('mu ~ dark*laser:');    disp(lme_mu.Coefficients);
disp('log(sigma^2) ~ dark*laser:'); disp(lme_sig.Coefficients);

% --- Direct subset comparisons: laser effect within each lighting level ----
% This is Option A from our discussion: filter to one lighting level, fit
% the single remaining binary factor. The resulting LaserOn coefficient
% *is* the laser effect in that lighting condition -- no beta arithmetic.
fprintf('\n  --- Laser effect within each lighting condition ---\n');
paramNames = {'Theta','Mu','LogSigma2'};
for lightLevel = ["dark","light"]
    subTbl = TrialParams(TrialParams.Dark == lightLevel, :);
    subTbl.LaserOn = removecats(subTbl.LaserOn);
    fprintf('  Laser effect in the %s (n=%d trials):\n', lightLevel, height(subTbl));
    for p = 1:numel(paramNames)
        lme_sub = fitlme(subTbl, sprintf('%s ~ LaserOn + (1|AnimalID)', paramNames{p}));
        rowIdx = strcmp(lme_sub.Coefficients.Name, 'LaserOn_laser_on');
        c = lme_sub.Coefficients(rowIdx, :);
        fprintf('    %-10s beta = %+.4f, SE = %.4f, p = %.4g\n', ...
            paramNames{p}, c.Estimate, c.SE, c.pValue);
    end
end

% --- Direct subset comparisons: lighting effect within each laser level ----
% Converse of the block above: filter to one laser level, fit the single
% remaining binary factor (Dark). The resulting Dark coefficient *is* the
% dark-vs-light effect in that laser condition -- no beta arithmetic.
fprintf('\n  --- Lighting effect within each laser condition ---\n');
paramNames = {'Theta','Mu','LogSigma2'};
for laserLevel = ["laser_off","laser_on"]
    subTbl = TrialParams(TrialParams.LaserOn == laserLevel, :);
    subTbl.Dark = removecats(subTbl.Dark);
    fprintf('  Dark vs light effect during %s (n=%d trials):\n', laserLevel, height(subTbl));
    for p = 1:numel(paramNames)
        lme_sub = fitlme(subTbl, sprintf('%s ~ Dark + (1|AnimalID)', paramNames{p}));
        rowIdx = strcmp(lme_sub.Coefficients.Name, 'Dark_dark');
        c = lme_sub.Coefficients(rowIdx, :);
        fprintf('    %-10s beta = %+.4f, SE = %.4f, p = %.4g\n', ...
            paramNames{p}, c.Estimate, c.SE, c.pValue);
    end
end

%% =========================================================================
%  PART B -- EVENT-TRIGGERED ANALYSIS (jump-diffusion / point-process link)
%% =========================================================================

%% B1. Peri-event traces (PSTH-style, for continuous z) --------------------
fprintf('\n=== [3/5] Part B1: peri-event traces ===\n');
nPre  = round(2.0*fs);   % 1 s before event
nPost = round(2.0*fs);   % 1 s after event
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
    plot(tAxis, mtrace, 'DisplayName', sprintf('%s (n=%d)', et, size(M,1)), 'linewidth', 3);
    fill([tAxis fliplr(tAxis)], [mtrace+setrace fliplr(mtrace-setrace)], 'k', ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility','off'        );
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
fprintf('  Fitting mixed model DeltaZ ~ EventType + (1|AnimalID) ...\n');
EventEffects.EventType = categorical(EventEffects.EventType);
lme_event = fitlme(EventEffects, 'DeltaZ ~ EventType + (1|AnimalID)');
disp(lme_event.Coefficients);
disp(anova(lme_event));

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


