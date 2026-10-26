%% =========================================================================
%  ANALYZE_ENGAGED_BOUT_SURVIVAL
%  Tests the "persistence" claim in the ZIm paper: does ZIm activation
%  lengthen ENGAGED bouts (follow + chase + hot pursuit, lumped), and is the
%  effect larger in the dark? Also checks whether the laser changes how
%  often engaged bouts START (re-engagement), to separate persistence from
%  initiation.
%
%  Differences from analyze_dwell_times.m (July):
%   - engaged states are lumped into bouts; moving between follow/chase/
%     hot pursuit continues a bout
%   - censored bouts are KEPT and handled by survival methods, not dropped
%   - bouts that end at capture / end of the laser window are CENSORED,
%     not counted as disengagements (cricket_present -> 0 at capture made
%     these look like ordinary bout ends before)
%   - inference respects n = 7 mice: Cox model stratified by mouse,
%     trial-label permutation test, and mouse-level bootstrap
%
%  INPUTS (workspace, Nframes x 1, as produced by the state detection +
%  overlap-resolution block in plotTPM4.m, i.e. AFTER overlaps are removed):
%     hotpursuit, chase, follow, stalk, wander, pause   (logical/0-1)
%     cricket_present, localframe, filename, metadata
%  metadata needs columns: filename, cricketdrop, laser_value, condition
%  (condition: 1 = light, 2 = dark, as in get_lightdark.m)
%
%  OUTPUTS: figure (6 panels), results struct, and two CSVs (bout table and
%  trial table) so the same bouts can be re-analyzed independently in Python.
% =========================================================================

%% 0. PARAMETERS -----------------------------------------------------------
fs              = 200;          % frames/s
windowSec       = 115;          % analysis window after cricket drop (s).
                                % Laser runs for 120 s starting 5 s BEFORE
                                % drop, so 115 s keeps laser trials fully
                                % stimulated; applied identically to
                                % laser-off trials.
gapTolSecList   = [0 0.1 0.25]; % bridge unclassified (NOTA) gaps <= this
                                % inside an engaged bout (sensitivity set)
primaryGapSec   = 0.1;          % primary analysis
tauSec          = 1.5;         % RMST horizon (s): longest at which every mouse x condition still has bouts at risk
nPerm           = 1000;          % trial-label permutations (increase for final)
nBoot           = 500;          % mouse-level bootstrap resamples
doPlot          = true;
rng(1);

statesOK = exist('hotpursuit','var') && exist('chase','var') && exist('follow','var') && ...
    exist('stalk','var') && exist('wander','var') && exist('pause','var') && ...
    exist('cricket_present','var') && exist('localframe','var') && ...
    exist('filename','var') && exist('metadata','var');
if ~statesOK
    error(['Missing inputs. Run the state detection and overlap-resolution ' ...
        'block (as in plotTPM4.m) first.']);
end
if ~exist('outputrootdir','var') || ~exist(outputrootdir,'dir'), outputrootdir = pwd; end
written = {};                   % every file this script writes, listed at the end

%% 1. TRIAL-LEVEL INFO (laser and light/dark are trial-level conditions) ----
fprintf('\n=== Trial info ===\n');
fnames = cellstr(filename(:));
[uf, ~, tid] = unique(fnames, 'stable');
tid = tid(:);
nT = numel(uf);
mdFile = cellstr(metadata.filename);

T_drop = nan(nT,1); T_laser = nan(nT,1); T_dark = nan(nT,1); T_animal = nan(nT,1);
for k = 1:nT
    r = find(cellfun(@(m) ~isempty(strfind(uf{k}, m)), mdFile), 1);
    if isempty(r), continue, end
    T_drop(k)  = metadata.cricketdrop(r);
    lv = metadata.laser_value(r); if isnan(lv), lv = 0; end
    T_laser(k) = double(lv ~= 0);
    T_dark(k)  = double(metadata.condition(r) == 2);
    tok = regexp(uf{k}, 'mouse-(\d+)', 'tokens', 'once');
    if ~isempty(tok), T_animal(k) = str2double(tok{1}); end
end
T_valid = ~isnan(T_drop) & T_drop > 0 & ~isnan(T_animal) & ~isnan(T_dark);
fprintf('  %d trials (%d valid), %d mice\n', nT, nnz(T_valid), numel(unique(T_animal(T_valid))));

%% 2. FRAME LABELS: 1 = engaged, -1 = disengaged, 0 = unclassified (NOTA),
%     NaN = outside analysis window (before drop, after capture or window end)
eng = logical(hotpursuit(:)) | logical(chase(:)) | logical(follow(:));
dis = (logical(stalk(:)) | logical(wander(:)) | logical(pause(:))) & ~eng;
tsd = (localframe(:) - T_drop(tid)) / fs;                     % time since drop
inWin = logical(cricket_present(:)) & tsd >= 0 & tsd <= windowSec & T_valid(tid);
L = nan(numel(eng),1);
L(inWin) = 0; L(inWin & dis) = -1; L(inWin & eng) = 1;
fprintf('  in-window frames: %.1f%% engaged, %.1f%% disengaged, %.1f%% unclassified\n', ...
    100*mean(L(inWin)==1), 100*mean(L(inWin)==-1), 100*mean(L(inWin)==0));

%% 3. EXTRACT BOUTS FOR EACH GAP TOLERANCE, FIT MODELS ---------------------
exitNames = {'direct disengage','NOTA -> disengage','NOTA -> re-engage','capture / window end'};
defNames  = {'primary (next labeled state disengaged)','any exit from engaged'};
R = struct();
for g = 1:numel(gapTolSecList)
    gapSec = gapTolSecList(g);
    B = local_extract_bouts(L, tid, round(gapSec*fs), fs, nT);
    B.animal = T_animal(B.trial); B.dark = T_dark(B.trial); B.laser = T_laser(B.trial);
    keep = ~B.leftTrunc;
    fn = fieldnames(B);
    for f = 1:numel(fn)
        if numel(B.(fn{f})) == numel(keep), B.(fn{f}) = B.(fn{f})(keep); end
    end
    evDef = {ismember(B.exitType, [1 2]), ismember(B.exitType, [1 2 3])};

    for d = 1:2
        for lt = 0:1                                   % lt = 0 light, 1 dark
            sel = B.dark == lt;
            [b, se, p] = local_cox(B.laser(sel), B.dur(sel), evDef{d}(sel), B.animal(sel));
            R(g).def(d).HR(lt+1)   = exp(b);
            R(g).def(d).CIlo(lt+1) = exp(b - 1.96*se);
            R(g).def(d).CIhi(lt+1) = exp(b + 1.96*se);
            R(g).def(d).pWald(lt+1) = p;
            R(g).def(d).beta(lt+1) = b;
        end
        % interaction model: laser, dark, laser x dark, stratified by mouse
        [bI, seI, pI] = local_cox([B.laser, B.dark, B.laser.*B.dark], B.dur, evDef{d}, B.animal);
        R(g).def(d).interactionBeta = bI(3);
        R(g).def(d).interactionP = pI(3);
    end
    R(g).gapSec = gapSec;
    R(g).B = B;
    R(g).evDef = evDef;
    fprintf('\n--- gap tolerance %.2f s: %d bouts ---\n', gapSec, numel(B.dur));
    for d = 1:2
        fprintf('  %s\n', defNames{d});
        fprintf('    laser HR, light: %.2f [%.2f %.2f], Wald p=%.3g\n', R(g).def(d).HR(1), R(g).def(d).CIlo(1), R(g).def(d).CIhi(1), R(g).def(d).pWald(1));
        fprintf('    laser HR, dark:  %.2f [%.2f %.2f], Wald p=%.3g\n', R(g).def(d).HR(2), R(g).def(d).CIlo(2), R(g).def(d).CIhi(2), R(g).def(d).pWald(2));
        fprintf('    laser x dark interaction p=%.3g\n', R(g).def(d).interactionP);
    end
end
gP = find(abs(gapTolSecList - primaryGapSec) < 1e-9, 1);
B = R(gP).B; ev = R(gP).evDef{1};

%% 4. EXIT COMPOSITION (diagnostic: if the laser changes HOW bouts end,
%     censoring choices could bias the hazard ratio) ----------------------
condIdx = 1 + B.dark*2 + B.laser;                     % 1 LOff 2 LOn 3 DOff 4 DOn
condNames = {'light off','light on','dark off','dark on'};
exitComp = zeros(4,4);
for c = 1:4
    for e = 1:4
        exitComp(c,e) = mean(B.exitType(condIdx==c) == e);
    end
end
fprintf('\n=== Exit composition (primary gap tolerance), rows = condition ===\n');
fprintf('%12s', ''); fprintf('%22s', exitNames{:}); fprintf('\n');
for c = 1:4, fprintf('%12s', condNames{c}); fprintf('%22.3f', exitComp(c,:)); fprintf('\n'); end

%% 5. INFERENCE FOR THE PRIMARY ANALYSIS -----------------------------------
% (a) permutation: laser was randomized at the trial level, so shuffle trial
%     laser labels within each mouse x lighting stratum (exact null of no
%     laser effect); refit stratified Cox in each lighting condition
fprintf('\n=== Permutation test (%d perms) ===\n', nPerm);
obsB = R(gP).def(1).beta;                             % [light dark]
permB = nan(nPerm, 2);
strat = T_animal*10 + T_dark;                         % trial-level strata
validTrials = find(T_valid);
for ip = 1:nPerm
    Lp = T_laser;
    us = unique(strat(validTrials));
    for s = us'
        idx = validTrials(strat(validTrials) == s);
        Lp(idx) = T_laser(idx(randperm(numel(idx))));
    end
    lasP = Lp(B.trial);
    for lt = 0:1
        sel = B.dark == lt;
        permB(ip, lt+1) = local_cox(lasP(sel), B.dur(sel), ev(sel), B.animal(sel));
    end
end
pPerm = (1 + sum(abs(permB) >= abs(obsB), 1)) ./ (nPerm + 1);
fprintf('  laser effect, light: HR %.2f, permutation p = %.3g\n', exp(obsB(1)), pPerm(1));
fprintf('  laser effect, dark:  HR %.2f, permutation p = %.3g\n', exp(obsB(2)), pPerm(2));

% (b) mouse-level bootstrap for CIs, including the dark - light difference
fprintf('\n=== Mouse bootstrap (%d resamples) ===\n', nBoot);
mice = unique(B.animal);
bootB = nan(nBoot, 2);
for ib = 1:nBoot
    pick = mice(randi(numel(mice), numel(mice), 1));
    idx = []; newStrat = [];
    for m = 1:numel(pick)
        ii = find(B.animal == pick(m));
        idx = [idx; ii]; newStrat = [newStrat; m*ones(numel(ii),1)]; %#ok<AGROW>
    end
    for lt = 0:1
        sel = B.dark(idx) == lt;
        if numel(unique(B.laser(idx(sel)))) < 2, continue, end
        bootB(ib, lt+1) = local_cox(B.laser(idx(sel)), B.dur(idx(sel)), ev(idx(sel)), newStrat(sel));
    end
end
bootHR = exp(bootB);
ciHR = prctile(bootHR, [2.5 97.5]);
diffB = bootB(:,2) - bootB(:,1);
ciDiff = prctile(exp(diffB), [2.5 97.5]);
fprintf('  laser HR light: %.2f, bootstrap 95%% CI [%.2f %.2f]\n', exp(obsB(1)), ciHR(1,1), ciHR(2,1));
fprintf('  laser HR dark:  %.2f, bootstrap 95%% CI [%.2f %.2f]\n', exp(obsB(2)), ciHR(1,2), ciHR(2,2));
fprintf('  ratio of HRs (dark/light): %.2f, bootstrap 95%% CI [%.2f %.2f]\n', exp(obsB(2)-obsB(1)), ciDiff(1), ciDiff(2));

%% 6. PER-MOUSE RESTRICTED MEAN SURVIVAL TIME (RMST, up to tauSec) ---------
fprintf('\n=== Per-mouse RMST (tau = %.1f s); bout duration quantiles: ', tauSec);
fprintf('50%% %.2f s, 90%% %.2f s, 99%% %.2f s ===\n', prctile(B.dur, [50 90 99]));
rmst = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        sel = B.animal == mice(m) & condIdx == c;
        if nnz(sel) >= 5, rmst(m,c) = local_rmst(B.dur(sel), ev(sel), tauSec); end
    end
end
[~, pRmstL] = local_signrank(rmst(:,2) - rmst(:,1));
[~, pRmstD] = local_signrank(rmst(:,4) - rmst(:,3));
fprintf('  light: off %.3f s, on %.3f s (paired signed-rank p = %.3g, n = %d mice)\n', ...
    mean(rmst(:,1),'omitnan'), mean(rmst(:,2),'omitnan'), pRmstL, nnz(all(~isnan(rmst(:,1:2)),2)));
fprintf('  dark:  off %.3f s, on %.3f s (paired signed-rank p = %.3g, n = %d mice)\n', ...
    mean(rmst(:,3),'omitnan'), mean(rmst(:,4),'omitnan'), pRmstD, nnz(all(~isnan(rmst(:,3:4)),2)));

%% 7. ONSET RATE: engaged bouts started per second of non-engaged time -----
onCount = accumarray(B.trial, 1, [nT 1]);
nonEngSec = accumarray(tid(L == 0 | L == -1), 1, [nT 1]) / fs;
useT = T_valid & nonEngSec > 1;
onRate = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        cT = 1 + T_dark*2 + T_laser;
        sel = useT & T_animal == mice(m) & cT == c;
        if any(sel), onRate(m,c) = sum(onCount(sel)) / sum(nonEngSec(sel)); end
    end
end
RRon = nan(1,2); pOn = nan(1,2);
if exist('fitglme', 'file') && exist('table', 'file')
    TT = table(onCount(useT), log(nonEngSec(useT)), T_laser(useT), T_dark(useT), ...
        categorical(T_animal(useT)), 'VariableNames', {'n','logT','Laser','Dark','Animal'});
    for lt = 0:1
        sub = TT(TT.Dark == lt, :);
        gm = fitglme(sub, 'n ~ Laser + (1|Animal)', 'Distribution', 'Poisson', 'Offset', sub.logT);
        row = strcmp(gm.Coefficients.Name, 'Laser');
        RRon(lt+1) = exp(gm.Coefficients.Estimate(row)); pOn(lt+1) = gm.Coefficients.pValue(row);
    end
    fprintf('\n=== Onset rate (Poisson GLMM, offset = non-engaged time) ===\n');
    fprintf('  laser rate ratio, light: %.2f (p = %.3g)\n', RRon(1), pOn(1));
    fprintf('  laser rate ratio, dark:  %.2f (p = %.3g)\n', RRon(2), pOn(2));
end

%% 8. FIGURE ---------------------------------------------------------------
if doPlot
    cOff = [0 0 0]; cOn = [0 0.6 0.9];
    fig = figure('Name', 'Engaged bout survival', 'Color', 'w', 'Position', [100 100 1800 1050]);
    set(fig, 'DefaultAxesFontSize', 18);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    lightNames = {'Light','Dark'};

    % A, B: Kaplan-Meier survival of engaged bouts
    for lt = 0:1
        nexttile; hold on;
        for la = 0:1
            sel = B.dark == lt & B.laser == la;
            [tt, S] = local_km(B.dur(sel), ev(sel));
            col = cOff; if la, col = cOn; end
            stairs([0; tt], [1; S], 'Color', col, 'LineWidth', 2.5);
        end
        xlim([0 tauSec]); ylim([0 1]); box off;
        xlabel('Engaged bout duration (s)'); ylabel('P(still engaged)');
        title(sprintf('%s: HR %.2f, perm p = %.2g', lightNames{lt+1}, exp(obsB(lt+1)), pPerm(lt+1)));
        legend({sprintf('laser off (%d)', nnz(B.dark==lt & B.laser==0)), ...
                sprintf('laser on (%d)',  nnz(B.dark==lt & B.laser==1))}, 'Location', 'northeast', 'Box', 'off');
    end

    % C: forest plot of laser hazard ratios across gap tolerances and event definitions
    nexttile; hold on;
    y = 0; yv = []; yl = {};
    for lt = 0:1
        for d = 1:2
            for g = 1:numel(gapTolSecList)
                y = y + 1;
                col = cOff; if lt, col = cOn; end
                mk = 'o'; if d == 2, mk = 's'; end
                plot([R(g).def(d).CIlo(lt+1) R(g).def(d).CIhi(lt+1)], [y y], '-', 'Color', col, 'LineWidth', 2);
                plot(R(g).def(d).HR(lt+1), y, mk, 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 9);
                tag = ''; if g == gP && d == 1, tag = ' *'; end
                yv(end+1) = y; %#ok<SAGROW>
                yl{end+1} = sprintf('%s, %s, gap %.2f%s', lightNames{lt+1}, char('A'+d-1), gapTolSecList(g), tag); %#ok<SAGROW>
            end
        end
        y = y + 1;   % blank row between light and dark
    end
    xline(1, '--k');
    set(gca, 'YTick', yv, 'YTickLabel', yl, 'YDir', 'reverse', 'FontSize', 18);
    set(gca, 'XScale', 'log'); xlabel('Laser hazard ratio (<1 = longer bouts)', 'FontSize', 18);
    title({'Sensitivity (Wald 95% CI)', 'A = primary, B = any exit; * = primary'}, 'FontSize', 18);
    box off;

    % D: per-mouse RMST
    nexttile; hold on;
    xs = [1 2; 4 5];
    for lt = 0:1
        cc = lt*2 + (1:2);
        plot(xs(lt+1,:), rmst(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(xs(lt+1,1), rmst(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
        plot(xs(lt+1,2), rmst(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
    end
    set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]);
    ylabel(sprintf('RMST to %.0f s (s)', tauSec));
    title(sprintf('Per mouse: light p=%.2g, dark p=%.2g', pRmstL, pRmstD)); box off;

    % E: per-mouse onset rate
    nexttile; hold on;
    for lt = 0:1
        cc = lt*2 + (1:2);
        plot(xs(lt+1,:), onRate(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(xs(lt+1,1), onRate(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
        plot(xs(lt+1,2), onRate(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
    end
    set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]);
    ylabel('Bout onsets / s non-engaged');
    title(sprintf('Onset rate: RR light %.2f, dark %.2f', RRon(1), RRon(2))); box off;

    % F: exit composition
    nexttile;
    bh = bar(exitComp, 'stacked');
    exitCols = [0.80 0.20 0.20; 0.95 0.60 0.50; 0.60 0.60 0.60; 0.30 0.60 0.30];
    for e = 1:4, bh(e).FaceColor = exitCols(e,:); end
    set(gca, 'XTickLabel', condNames); ylabel('Fraction of bouts'); ylim([0 1]);
    legend(exitNames, 'Location', 'southoutside', 'Box', 'off', 'FontSize', 18);
    title('How engaged bouts end'); box off;

    title(tl, sprintf('Engaged bouts (follow + chase + hot pursuit), gap tolerance %.2f s, first %d s after drop', ...
        primaryGapSec, windowSec), 'FontSize', 18);
    outFile = fullfile(outputrootdir, 'engaged_bout_survival.pdf'); exportgraphics(fig, outFile, 'ContentType', 'vector'); written{end+1} = outFile;
end

%% 9. SAVE RESULTS AND EXPORT BOUT / TRIAL TABLES FOR INDEPENDENT CHECKS ----
results = struct('R', R, 'pPerm', pPerm, 'bootHR', bootHR, 'ciHR', ciHR, 'ciDiff', ciDiff, ...
    'rmst', rmst, 'onRate', onRate, 'RRon', RRon, 'pOn', pOn, 'exitComp', exitComp, ...
    'mice', mice, 'params', struct('fs', fs, 'windowSec', windowSec, 'gapTolSecList', gapTolSecList, ...
    'primaryGapSec', primaryGapSec, 'tauSec', tauSec, 'nPerm', nPerm, 'nBoot', nBoot));
outFile = fullfile(outputrootdir, 'engaged_bout_survival_results.mat'); save(outFile, 'results'); written{end+1} = outFile;
if exist('writetable', 'file') && exist('table', 'file')
    BT = table(B.trial, uf(B.trial), B.animal, B.dark, B.laser, B.startFrame, B.dur, B.exitType, ev, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','start_frame','dur_s','exit_type','event_primary'});
    outFile = fullfile(outputrootdir, 'engaged_bouts_primary.csv'); writetable(BT, outFile); written{end+1} = outFile;
    TTr = table((1:nT)', uf, T_animal, T_dark, T_laser, T_drop, onCount, nonEngSec, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','cricketdrop','n_bout_onsets','nonengaged_s'});
    outFile = fullfile(outputrootdir, 'engaged_bout_trials_primary.csv'); writetable(TTr, outFile); written{end+1} = outFile;
end
fprintf('\nFiles written:\n'); fprintf('  %s\n', written{:}); fprintf('Done.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function B = local_extract_bouts(L, tid, gapFrames, fs, nT)
% Run-length encode the label vector within trials, bridge short NOTA gaps
% between engaged runs, and classify how each engaged bout ends.
% exitType: 1 = next frame disengaged; 2 = NOTA, then disengaged;
%           3 = NOTA (longer than gap tolerance), then engaged again;
%           4 = capture / end of analysis window / end of trial (censored)
v = L; v(isnan(v)) = -9;
[rs, re, rv, rt] = local_runs(v, tid);
if gapFrames > 0
    len = re - rs + 1;
    in = (2:numel(rs)-1)';
    br = in(rv(in) == 0 & len(in) <= gapFrames & rv(in-1) == 1 & rv(in+1) == 1 & ...
        rt(in-1) == rt(in) & rt(in+1) == rt(in));
    for i = br', v(rs(i):re(i)) = 1; end
    [rs, re, rv, rt] = local_runs(v, tid);
end
nr = numel(rs);
ib = find(rv == 1);
exitType = 4*ones(numel(ib), 1);
leftTrunc = false(numel(ib), 1);
for j = 1:numel(ib)
    i = ib(j);
    leftTrunc(j) = (i == 1) || rt(i-1) ~= rt(i) || rv(i-1) == -9;  % engaged at window start
    if i == nr || rt(i+1) ~= rt(i) || rv(i+1) == -9
        exitType(j) = 4;
    elseif rv(i+1) == -1
        exitType(j) = 1;
    elseif i+2 <= nr && rt(i+2) == rt(i) && rv(i+2) == -1
        exitType(j) = 2;
    elseif i+2 <= nr && rt(i+2) == rt(i) && rv(i+2) == 1
        exitType(j) = 3;
    else
        exitType(j) = 4;
    end
end
B.trial = rt(ib);
B.startFrame = rs(ib);
B.dur = (re(ib) - rs(ib) + 1) / fs;
B.exitType = exitType;
B.leftTrunc = leftTrunc;
B.nT = nT;
end

function [rs, re, rv, rt] = local_runs(v, tid)
n = numel(v);
chg = [true; diff(v) ~= 0 | diff(tid) ~= 0];
rs = find(chg); re = [rs(2:end) - 1; n];
rv = v(rs); rt = tid(rs);
end

function [b, se, p] = local_cox(X, t, ev, strata)
% Stratified Cox proportional hazards (Efron ties); censoring = ~event
[b, ~, ~, st] = coxphfit(X, t, 'Censoring', ~logical(ev), 'Strata', strata, 'Ties', 'efron');
se = st.se; p = st.p;
end

function [tt, S] = local_km(t, ev)
% Kaplan-Meier survivor function; tt = distinct event times, S = S(tt)
t = t(:); ev = logical(ev(:));
tt = unique(t(ev));
if isempty(tt), S = []; return, end
nAt = sum(t >= tt', 1)';
dAt = sum((t == tt') & ev, 1)';
S = cumprod(1 - dAt ./ nAt);
end

function r = local_rmst(t, ev, tau)
% Restricted mean survival time: area under the KM curve from 0 to tau
[tt, S] = local_km(t, ev);
k = tt <= tau; tt = tt(k); S = S(k);
r = sum([1; S] .* diff([0; tt; tau]));
end

function [h, p] = local_signrank(d)
d = d(~isnan(d));
if numel(d) < 2, h = false; p = NaN; return, end
p = signrank(d); h = p < 0.05;
end
