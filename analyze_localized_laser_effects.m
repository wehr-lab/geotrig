%% =========================================================================
%  ANALYZE_LOCALIZED_LASER_EFFECTS
%  Generic battery of tests for a spatially localized (not simply mean- or
%  variance-shifting) effect of a binary factor on the latent engagement
%  axis z(t). Runs the battery for both:
%    (1) laser-on vs laser-off, within dark trials
%    (2) dark vs light, within laser-off trials
%  then formally tests the "mirror hypothesis" -- that the laser-on effect
%  (in the dark) is essentially the same as the light-vs-dark effect (with
%  laser off) -- by comparing per-bin drift effects on a SHARED z-axis grid.
%
%  PART 1 (per run): trial-level occupancy at each extreme (binomial GLMM)
%  PART 2 (per run): frame-level binned drift/diffusion (Kramers-Moyal)
%  PART 3: cross-comparison of the two runs' per-bin drift effects, on
%          shared bin edges, to test the mirror hypothesis directly.
%  PART 4: checks where each run's peak-effect bin actually sits, against
%          ALL SIX behavioral states' z-ranges (not just the one state
%          that happened to show up in the per-state corroborating check --
%          that check can be misled by spillover from a broad effect).
% =========================================================================

%% 0. ASSUMED INPUTS (frame-level, all Nx1, same length) ---------------------
% z, trial_id, animal_id, state_id, dark, laser_on : as used throughout
% fs : frame rate
fs=fps;
laser_on=laseron;

set(groot, 'defaultAxesFontSize', 18);

%% Shared bin edges, computed ONCE from all valid frames -----------------
% Using one shared grid (rather than each run picking its own bins from its
% own subset) is what makes "bin 6" mean the same region of z in both runs,
% so their per-bin effects can be directly compared below.
nBins = 10;
zAll = z(~isnan(z));
sharedEdges = prctile(zAll, linspace(0, 100, nBins+1));
sharedEdges(1) = -inf; sharedEdges(end) = inf;
if numel(unique(sharedEdges)) < numel(sharedEdges)
    warning('Duplicate shared bin edges (ties in z) -- reducing bin count.');
    sharedEdges = unique(sharedEdges);
end

%% Run 1: laser effect, within dark trials -----------------------------------
laserResults = local_runLocalizedEffectAnalysis(z, trial_id, animal_id, state_id, fs, ...
    dark, 'dark trials', laser_on, {'laser_off','laser_on'}, 'LaserOn', sharedEdges);

%% Run 2: lighting effect, within laser-off trials ---------------------------
darkResults = local_runLocalizedEffectAnalysis(z, trial_id, animal_id, state_id, fs, ...
    ~laser_on, 'laser-off trials', dark, {'light','dark'}, 'Dark', sharedEdges);

%% =========================================================================
%  PART 3 -- TESTING THE MIRROR HYPOTHESIS
%  H: LaserOn effect (in dark) ~= Light-vs-Dark effect (laser off)
%     = -1 * (Dark-vs-Light effect), since "dark" was the printed contrast
% =========================================================================
fprintf('\n########################################################\n');
fprintf('### Testing mirror hypothesis: laser-on (dark) vs light-vs-dark (laser off)\n');
fprintf('########################################################\n');

laserBeta = laserResults.beta;      % effect of laser-on, per bin, in dark
laserSE   = laserResults.SE;
lightBeta = -darkResults.beta;      % NEGATE: darkResults.beta is dark-vs-light;
lightSE   = darkResults.SE;         % light-vs-dark is its negation. SE unaffected by sign flip.
binCenters = laserResults.binCenters;   % same grid in both runs by construction

fprintf('\n  %-6s %-8s %-14s %-14s %-14s\n', 'bin', 'z', 'laser beta', 'light beta', 'diff (z-test)');
diffZ = nan(size(binCenters)); diffP = nan(size(binCenters));
for b = 1:numel(binCenters)
    if isnan(laserBeta(b)) || isnan(lightBeta(b)), continue; end
    d = laserBeta(b) - lightBeta(b);
    seD = sqrt(laserSE(b)^2 + lightSE(b)^2);
    diffZ(b) = d / seD;
    diffP(b) = 2*(1 - normcdf(abs(diffZ(b))));
    fprintf('  %-6d %-8.2f %-14.4f %-14.4f z=%.2f, p=%.4g\n', ...
        b, binCenters(b), laserBeta(b), lightBeta(b), diffZ(b), diffP(b));
end

validBoth = ~isnan(laserBeta) & ~isnan(lightBeta);
[rho, pRho] = corr(laserBeta(validBoth), lightBeta(validBoth), 'Type', 'Pearson');
fprintf('\n  Cross-bin correlation (laser beta vs light beta): r = %.3f, p = %.4g (n = %d bins)\n', ...
    rho, pRho, nnz(validBoth));
fprintf('  If the mirror hypothesis holds: expect r close to +1 (same shape, same sign,\n');
fprintf('  similar magnitude) and few/no bins with a significant diff (z-test) above.\n');
fprintf('  A weak/non-significant r, or a cluster of significant per-bin differences,\n');
fprintf('  is evidence the two manipulations act on different parts of the range.\n');

figure('Name', 'Mirror hypothesis: laser vs light effect, per bin');
subplot(1,2,1); hold on;
plot(binCenters, laserBeta, '-o', 'DisplayName', 'laser-on effect (dark)');
plot(binCenters, lightBeta, '-o', 'DisplayName', 'light-vs-dark effect (laser off)');
yline(0, '--k', 'HandleVisibility','off');
xlabel('z (engagement axis)', 'FontSize', 18); ylabel('per-bin drift beta', 'FontSize', 18);
title('Per-bin effect comparison', 'FontSize', 18);
legend('Location','best', 'FontSize', 18); grid on;

subplot(1,2,2); hold on;
scatter(lightBeta, laserBeta, 60, 'filled');
lims = [min([laserBeta;lightBeta],[],'omitnan'), max([laserBeta;lightBeta],[],'omitnan')];
plot(lims, lims, '--k', 'DisplayName', 'perfect mirror (y=x)');
for b = 1:numel(binCenters)
    if validBoth(b)
        text(lightBeta(b), laserBeta(b), sprintf(' %d', b), 'FontSize', 18);
    end
end
xlabel('light-vs-dark beta', 'FontSize', 18); ylabel('laser-on beta', 'FontSize', 18);
title(sprintf('r = %.2f, p = %.3g', rho, pRho), 'FontSize', 18);
legend('Location','best', 'FontSize', 18); grid on; axis equal;

%% =========================================================================
%  PART 4 -- WHERE DO THE DIVERGENT BINS (6 AND 9) ACTUALLY SIT?
%  Checks bin 6 (laser's peak effect) and bin 9 (light's peak effect)
%  against ALL SIX behavioral states' z-ranges, not just chase -- an
%  earlier chase-only version of this check turned out to be misleading:
%  bin 9 does sit inside chase, but bin 6 sits well BELOW chase entirely,
%  meaning chase showing up as significant for the laser effect in the
%  per-state check was likely spillover from a broad effect rather than
%  evidence that laser's mechanism is centered in chase.
% =========================================================================
fprintf('\n########################################################\n');
fprintf('### Where do bins 6 and 9 actually sit, across all six states?\n');
fprintf('########################################################\n');

statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
bin6z = laserResults.binCenters(6);
bin9z = laserResults.binCenters(9);

stateLo = nan(6,1); stateHi = nan(6,1); stateMed = nan(6,1);
for s = 1:6
    sel = (state_id == s) & ~isnan(z);
    if nnz(sel) < 20
        fprintf('  %-12s insufficient data, skipped\n', statenames{s});
        continue
    end
    stateLo(s) = prctile(z(sel), 5);
    stateHi(s) = prctile(z(sel), 95);
    stateMed(s) = median(z(sel), 'omitnan');
    fprintf('  %-12s z-range (5-95pct): [%.3f, %.3f], median = %.3f -- bin6 %s, bin9 %s\n', ...
        statenames{s}, stateLo(s), stateHi(s), stateMed(s), ...
        local_insideStr(bin6z, stateLo(s), stateHi(s)), ...
        local_insideStr(bin9z, stateLo(s), stateHi(s)));
end

%% Density-based assignment -- sharper than range containment when several
%% states' ranges overlap (as stalk/wander/pause do here). At each bin
%% center, evaluate each state's kernel density and rank states by it,
%% rather than just checking whether the point falls inside a broad range.
fprintf('\n  --- Density-ranked state assignment (sharper than range containment) ---\n');
for whichBin = ["bin6","bin9"]
    if whichBin == "bin6", zPoint = bin6z; label = 'Bin 6 (laser peak)';
    else,                  zPoint = bin9z; label = 'Bin 9 (light peak)'; end

    dens = nan(6,1);
    for s = 1:6
        sel = (state_id == s) & ~isnan(z);
        if nnz(sel) < 20, continue; end
        dens(s) = ksdensity(z(sel), zPoint);
    end
    [sortedDens, order] = sort(dens, 'descend', 'MissingPlacement','last');
    fprintf('  %s, z=%+.3f:\n', label, zPoint);
    for r = 1:6
        s = order(r);
        if isnan(sortedDens(r)), continue; end
        fprintf('    %d. %-12s density = %.4f\n', r, statenames{s}, sortedDens(r));
    end
end

figure('Name', 'Bin 6 / bin 9 relative to all state z-ranges');
hold on;
colors = lines(6);
for s = 1:6
    sel = (state_id == s) & ~isnan(z);
    if nnz(sel) < 20, continue; end
    histogram(z(sel), 60, 'Normalization', 'pdf', 'FaceColor', colors(s,:), ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', statenames{s});
end
xline(bin6z, '-r', 'LineWidth', 3, 'DisplayName', 'bin 6 (laser peak)');
xline(bin9z, '-b', 'LineWidth', 3, 'DisplayName', 'bin 9 (light peak)');
xlabel('z (engagement axis)', 'FontSize', 18); ylabel('density', 'FontSize', 18);
title('All state z-distributions vs. bin 6 / bin 9 locations', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 18); grid on;

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function s = local_insideStr(v, lo, hi)
    if v >= lo && v <= hi
        s = sprintf('INSIDE [%.3f, %.3f]', lo, hi);
    else
        s = sprintf('outside [%.3f, %.3f]', lo, hi);
    end
end

function results = local_runLocalizedEffectAnalysis(z, trial_id, animal_id, state_id, fs, ...
    subsetMask, subsetLabel, compareVar, compareLabels, compareVarName, sharedEdges)
    % Runs the occupancy GLMM + binned drift/diffusion + per-state
    % corroborating check, testing whether `compareVar` (binary) has a
    % spatially localized effect on z(t), within frames where subsetMask
    % is true. Uses sharedEdges (same grid across runs) rather than
    % recomputing bin edges from this subset alone, so results are
    % directly comparable to another run of this function.
    %
    % compareLabels = {referenceLevelName, otherLevelName} -- the SECOND
    % label is the one whose coefficient gets extracted/printed.
    %
    % Returns: results.binCenters, results.beta, results.SE, results.pValue
    % (per bin, for the drift comparison), for cross-run comparison.

    dt = 1/fs;
    inSubset = logical(subsetMask);
    cmp = logical(compareVar);
    betaName = sprintf('CompareVar_%s', compareLabels{2});

    fprintf('\n########################################################\n');
    fprintf('### Localized effect of %s, within %s\n', compareVarName, subsetLabel);
    fprintf('########################################################\n');

    %% Part 1: trial-level occupancy at each extreme (binomial GLMM) -------
    fprintf('=== Part 1: occupancy at each extreme (%s only) ===\n', subsetLabel);

    zSub = z(inSubset & ~isnan(z));
    hiThresh = prctile(zSub, 90);
    loThresh = prctile(zSub, 10);
    fprintf('  Hot-pursuit-end threshold (90th pct): z > %.3f\n', hiThresh);
    fprintf('  Pause-end threshold (10th pct):       z < %.3f\n', loThresh);

    trial_ids = unique(trial_id(inSubset & ~isnan(trial_id)));
    nT = numel(trial_ids);
    TrialID = trial_ids(:); AnimalID = nan(nT,1); CompareVar = false(nT,1);
    NTotal = nan(nT,1); NAboveHi = nan(nT,1); NBelowLo = nan(nT,1);

    nbytes = fprintf('  %.1f%%', 0);
    for i = 1:nT
        idx = (trial_id == trial_ids(i)) & inSubset;
        zseg = z(idx);
        valid = ~isnan(zseg);
        NTotal(i) = nnz(valid);
        NAboveHi(i) = nnz(zseg(valid) > hiThresh);
        NBelowLo(i) = nnz(zseg(valid) < loThresh);

        cmpVals = cmp(idx); cmpVals = cmpVals(~isnan(double(cmpVals)));
        if ~isempty(cmpVals), CompareVar(i) = cmpVals(1); end
        animVals = animal_id(idx); animVals = animVals(~isnan(animVals));
        if ~isempty(animVals), AnimalID(i) = animVals(1); end

        fprintf(repmat('\b',1,nbytes));
        nbytes = fprintf('  %.1f%%', 100*i/nT);
    end
    fprintf('\n');

    OccTbl = table(TrialID, AnimalID, CompareVar, NTotal, NAboveHi, NBelowLo);
    OccTbl(OccTbl.NTotal < 30, :) = [];
    OccTbl.CompareVar = categorical(OccTbl.CompareVar, [false true], compareLabels);

    fprintf('\n  --- Hot-pursuit-end occupancy ~ %s (n=%d trials) ---\n', compareVarName, height(OccTbl));
    glme_hi = fitglme(OccTbl, 'NAboveHi ~ CompareVar + (1|AnimalID)', ...
        'Distribution', 'binomial', 'BinomialSize', OccTbl.NTotal);
    disp(glme_hi.Coefficients);

    fprintf('\n  --- Pause-end occupancy ~ %s (n=%d trials) ---\n', compareVarName, height(OccTbl));
    glme_lo = fitglme(OccTbl, 'NBelowLo ~ CompareVar + (1|AnimalID)', ...
        'Distribution', 'binomial', 'BinomialSize', OccTbl.NTotal);
    disp(glme_lo.Coefficients);

    %% Part 2: binned drift/diffusion (Kramers-Moyal), on SHARED bins ------
    fprintf('\n=== Part 2: binned drift/diffusion along z (%s only) ===\n', subsetLabel);

    isNextSameTrial = trial_id(1:end-1) == trial_id(2:end);
    validPair = isNextSameTrial & ~isnan(z(1:end-1)) & ~isnan(z(2:end)) ...
                & inSubset(1:end-1) & inSubset(2:end);
    idxValid = find(validPair);

    zCur         = z(idxValid);
    dz           = z(idxValid + 1) - z(idxValid);
    cmpAtStep    = logical(cmp(idxValid));
    animalAtStep = animal_id(idxValid);

    edges = sharedEdges;
    nBins = numel(edges) - 1;
    binIdx = discretize(zCur, edges);

    drift = nan(nBins,2); driftSE = nan(nBins,2);
    diffusion = nan(nBins,2); diffusionSE = nan(nBins,2);
    binCenters = nan(nBins,1);

    nbytes = fprintf('  %.1f%%', 0);
    step = 0; totalSteps = nBins*2;
    for b = 1:nBins
        binCenters(b) = median(zCur(binIdx==b), 'omitnan');
        for c = 1:2
            isOther = (c==2);
            sel = (binIdx==b) & (cmpAtStep==isOther);
            n = nnz(sel);
            if n >= 20
                incrRate = dz(sel) / dt;
                drift(b,c) = mean(incrRate);
                driftSE(b,c) = std(incrRate) / sqrt(n);
                diffusion(b,c) = var(dz(sel)) / dt;
                diffusionSE(b,c) = diffusion(b,c) * sqrt(2/(n-1));
            end
            step = step + 1;
            fprintf(repmat('\b',1,nbytes));
            nbytes = fprintf('  %.1f%%', 100*step/totalSteps);
        end
    end
    fprintf('\n');

    figure('Name', sprintf('%s effect within %s', compareVarName, subsetLabel));
    subplot(2,1,1); hold on;
    errorbar(binCenters, drift(:,1), driftSE(:,1), '-o', 'DisplayName', compareLabels{1});
    errorbar(binCenters, drift(:,2), driftSE(:,2), '-o', 'DisplayName', compareLabels{2});
    yline(0, '--k', 'HandleVisibility','off');
    local_overlayStateRanges(z, state_id, inSubset);
    xlabel('z (engagement axis)', 'FontSize', 18); ylabel('drift, E[\Deltaz]/\Deltat', 'FontSize', 18);
    title(sprintf('Binned drift along z -- %s (%s)', compareVarName, subsetLabel), 'FontSize', 18);
    legend('Location','best', 'FontSize', 18); grid on;

    subplot(2,1,2); hold on;
    errorbar(binCenters, diffusion(:,1), diffusionSE(:,1), '-o', 'DisplayName', compareLabels{1});
    errorbar(binCenters, diffusion(:,2), diffusionSE(:,2), '-o', 'DisplayName', compareLabels{2});
    local_overlayStateRanges(z, state_id, inSubset);
    xlabel('z (engagement axis)', 'FontSize', 18); ylabel('diffusion, Var[\Deltaz]/\Deltat', 'FontSize', 18);
    title(sprintf('Binned diffusion along z -- %s (%s)', compareVarName, subsetLabel), 'FontSize', 18);
    legend('Location','best', 'FontSize', 18); grid on;

    %% Per-bin significance, and collect for cross-run comparison ----------
    fprintf('\n  --- Per-bin drift comparison (%s) ---\n', compareVarName);
    beta = nan(nBins,1); SE = nan(nBins,1); pValue = nan(nBins,1);
    for b = 1:nBins
        sel = (binIdx == b);
        if nnz(sel & cmpAtStep) < 20 || nnz(sel & ~cmpAtStep) < 20
            fprintf('  bin %2d (z~%+.2f): insufficient data, skipped\n', b, binCenters(b));
            continue
        end
        BinTbl = table(dz(sel)/dt, ...
            categorical(cmpAtStep(sel), [false true], compareLabels), ...
            animalAtStep(sel), ...
            'VariableNames', {'IncrRate','CompareVar','AnimalID'});
        lme_bin = fitlme(BinTbl, 'IncrRate ~ CompareVar + (1|AnimalID)');
        rowIdx = strcmp(lme_bin.Coefficients.Name, betaName);
        c = lme_bin.Coefficients(rowIdx, :);
        beta(b) = c.Estimate; SE(b) = c.SE; pValue(b) = c.pValue;
        fprintf('  bin %2d (z~%+.2f): beta=%+.4f, SE=%.4f, p=%.4g\n', ...
            b, binCenters(b), beta(b), SE(b), pValue(b));
    end

    %% Per-bin diffusion (variance) comparison -- Brown-Forsythe style -------
    % Standard Levene/Brown-Forsythe test: compare |increment - group median|
    % between groups (median-based, not mean-based, for robustness to
    % skewed increment distributions), fit with a mixed model so animal
    % clustering is still accounted for, exactly like the drift test above.
    % A positive beta means the SECOND compareLabels level has MORE spread
    % (larger diffusion) than the reference level in that bin.
    fprintf('\n  --- Per-bin diffusion comparison (%s) ---\n', compareVarName);
    diffBeta = nan(nBins,1); diffSE = nan(nBins,1); diffPValue = nan(nBins,1);
    for b = 1:nBins
        sel = (binIdx == b);
        if nnz(sel & cmpAtStep) < 20 || nnz(sel & ~cmpAtStep) < 20
            fprintf('  bin %2d (z~%+.2f): insufficient data, skipped\n', b, binCenters(b));
            continue
        end
        grp = cmpAtStep(sel);
        dzSel = dz(sel);
        adDev = nan(size(dzSel));
        adDev(grp)  = abs(dzSel(grp)  - median(dzSel(grp),  'omitnan'));
        adDev(~grp) = abs(dzSel(~grp) - median(dzSel(~grp), 'omitnan'));

        DiffTbl = table(adDev, ...
            categorical(grp, [false true], compareLabels), ...
            animalAtStep(sel), ...
            'VariableNames', {'AbsDev','CompareVar','AnimalID'});
        lme_diff = fitlme(DiffTbl, 'AbsDev ~ CompareVar + (1|AnimalID)');
        rowIdx = strcmp(lme_diff.Coefficients.Name, betaName);
        c = lme_diff.Coefficients(rowIdx, :);
        diffBeta(b) = c.Estimate; diffSE(b) = c.SE; diffPValue(b) = c.pValue;
        fprintf('  bin %2d (z~%+.2f): beta=%+.4f, SE=%.4f, p=%.4g  (positive = %s has larger spread)\n', ...
            b, binCenters(b), diffBeta(b), diffSE(b), diffPValue(b), compareLabels{2});
    end

    %% Corroborating check: per-state drift comparison ----------------------
    fprintf('\n=== Corroborating check: drift by discrete behavioral state (%s) ===\n', compareVarName);
    statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
    stateAtStep = state_id(idxValid);

    for s = 1:numel(statenames)
        sel = (stateAtStep == s);
        if nnz(sel & cmpAtStep) < 20 || nnz(sel & ~cmpAtStep) < 20
            fprintf('  %-12s insufficient data, skipped\n', statenames{s});
            continue
        end
        StateTbl = table(dz(sel)/dt, ...
            categorical(cmpAtStep(sel), [false true], compareLabels), ...
            animalAtStep(sel), ...
            'VariableNames', {'IncrRate','CompareVar','AnimalID'});
        lme_state = fitlme(StateTbl, 'IncrRate ~ CompareVar + (1|AnimalID)');
        rowIdx = strcmp(lme_state.Coefficients.Name, betaName);
        c = lme_state.Coefficients(rowIdx, :);
        fprintf('  %-12s beta=%+.4f, SE=%.4f, p=%.4g (n=%d)\n', ...
            statenames{s}, c.Estimate, c.SE, c.pValue, nnz(sel));
    end

    fprintf('\nDone with %s within %s.\n', compareVarName, subsetLabel);

    results = struct('binCenters', binCenters, 'beta', beta, 'SE', SE, 'pValue', pValue, ...
        'diffBeta', diffBeta, 'diffSE', diffSE, 'diffPValue', diffPValue);
end

function local_overlayStateRanges(z, state_id, subsetMask)
    % Shades each behavioral state's 5th-95th percentile z-range as a
    % translucent vertical band, for visually checking whether effect
    % locations in the decile-bin analysis line up with a specific state.
    statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
    colors = lines(numel(statenames));
    yl = ylim;
    for s = 1:numel(statenames)
        sel = (state_id == s) & subsetMask & ~isnan(z);
        if nnz(sel) < 20, continue; end
        lo = prctile(z(sel), 5); hi = prctile(z(sel), 95);
        patch([lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], colors(s,:), ...
            'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        text(median(z(sel), 'omitnan'), yl(2), statenames{s}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 18, 'Color', colors(s,:)*0.7);
    end
    ylim(yl);
end
