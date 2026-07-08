%% =========================================================================
%  TEST_GEOMETRIC_VARIABLE_INDEPENDENCE
%  Trial-wise Pearson correlation between geometric variables (mouse_spd,
%  range, az), tested for a nonzero median via Wilcoxon signed-rank test.
%  This avoids the invalid p-value from a single pooled-frame correlation,
%  which treats thousands of autocorrelated frames as independent samples.
% =========================================================================

%% 0. ASSUMED INPUTS ----------------------------------------------------------
% mouse_spd, range, az, trial_id, animal_id : Nx1, frame-level
%
% NOTE on az: if az is a circular/signed bearing (wraps at +/-180 deg or
% +/-pi), Pearson correlation on the RAW variable can be misleading near
% the wrap point. If your Figure 1 az is already handled that way, this is
% fine; if not, consider testing against az_dev (|wrapped az|, as used in
% estimate_engagement_axis.m) instead, or a circular-correlation measure.

minFramesPerTrial = 30;   % minimum valid frames required to trust a trial's r

varNames = {'mouse_spd','range','az'};
varData  = {mouse_spd, range, az};
nVars = numel(varData);

trial_ids = unique(trial_id(~isnan(trial_id)));
nT = numel(trial_ids);

%% 1. Compute trial-level Pearson r for every variable pair ------------------
pairIdx = nchoosek(1:nVars, 2);   % all unique pairs, e.g. [1 2; 1 3; 2 3]
nPairs = size(pairIdx, 1);

trialR = nan(nT, nPairs);
trialAnimal = nan(nT, 1);

fprintf('Computing trial-level correlations for %d trials, %d variable pairs ...\n', nT, nPairs);
nbytes = fprintf('  %.1f%%', 0);
for i = 1:nT
    idx = (trial_id == trial_ids(i));
    for p = 1:nPairs
        x = varData{pairIdx(p,1)}(idx);
        y = varData{pairIdx(p,2)}(idx);
        valid = ~isnan(x) & ~isnan(y);
        if nnz(valid) >= minFramesPerTrial
            trialR(i,p) = corr(x(valid), y(valid), 'Type', 'Pearson');
        end
    end
    animVals = animal_id(idx); animVals = animVals(~isnan(animVals));
    if ~isempty(animVals), trialAnimal(i) = animVals(1); end

    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('  %.1f%%', 100*i/nT);
end
fprintf('\n');

%% 2. Wilcoxon signed-rank test per pair: is median trial-level r != 0? ------
fprintf('\n--- Wilcoxon signed-rank test: trial-level r vs. 0 ---\n');
fprintf('(primary test -- one r per trial, trial is the unit of replication)\n\n');
for p = 1:nPairs
    r = trialR(:,p); r = r(~isnan(r));
    medR = median(r);
    pval = signrank(r);
    fprintf('  %-10s vs %-10s: median r = %+.3f, n = %d trials, signrank p = %.4g\n', ...
        varNames{pairIdx(p,1)}, varNames{pairIdx(p,2)}, medR, numel(r), pval);
end
fprintf('\n  (3 pairs tested -- Bonferroni-corrected alpha would be 0.05/3 = 0.0167)\n');

%% 3. Robustness check: collapse to one median r per animal first -----------
% Guards against pseudoreplication if trial counts differ substantially
% across animals -- median trial-level r per animal, then Wilcoxon across
% animals instead of across trials. If this and the trial-level test above
% disagree, the trial-level result may be dominated by one or two
% heavily-sampled animals.
fprintf('\n--- Same test, collapsed to one median r per animal (robustness check) ---\n');
animals = unique(trialAnimal(~isnan(trialAnimal)));
for p = 1:nPairs
    animalR = nan(numel(animals), 1);
    for a = 1:numel(animals)
        sel = (trialAnimal == animals(a));
        animalR(a) = median(trialR(sel,p), 'omitnan');
    end
    animalR = animalR(~isnan(animalR));
    medR = median(animalR);
    pval = signrank(animalR);
    fprintf('  %-10s vs %-10s: median r = %+.3f, n = %d animals, signrank p = %.4g\n', ...
        varNames{pairIdx(p,1)}, varNames{pairIdx(p,2)}, medR, numel(animalR), pval);
end
