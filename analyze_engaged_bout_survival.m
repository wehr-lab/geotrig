%% =========================================================================
%  ANALYZE_ENGAGED_BOUT_SURVIVAL
%  Tests the "persistence" claim in the ZIm paper. PRIMARY measure is the
%  ENGAGEMENT EPISODE: a chain of engaged bouts (follow + chase + hot
%  pursuit, lumped) linked by unclassified interruptions that end in
%  renewed engagement. An episode ends when the mouse disengages (stalk,
%  wander or pause; event) or at capture / end of window (censored).
%
%  Questions, in order:
%   1. Does ZIm activation lengthen engagement episodes (persistence)?
%      Episode survival: KM, Cox stratified by mouse, trial-label
%      permutation, mouse bootstrap, per-mouse RMST.
%   2. Do continuous engaged bouts lengthen? (negative control)
%      Bout survival with ANY exit from engagement as the event.
%   3. Does the laser start more episodes (initiation)?
%      Episodes started per second outside episodes.
%   4. How do bouts end? (exit composition; explains 1 vs 2)
%
%  Exit types of an engaged bout:
%    1 direct disengage  2 NOTA -> disengage  3 NOTA -> re-engage
%    4 capture / window end (censored)
%  Episodes chain bouts across exit type 3.
%
%  INPUTS (workspace, Nframes x 1, as produced by the state detection +
%  overlap-resolution block in plotTPM4.m, i.e. AFTER overlaps are removed):
%     hotpursuit, chase, follow, stalk, wander, pause   (logical/0-1)
%     cricket_present, localframe, filename, metadata
%  metadata needs columns: filename, cricketdrop, laser_value, condition
%  (condition: 1 = light, 2 = dark, as in get_lightdark.m)
%
%  OUTPUTS: figure (8 panels), results struct, and three CSVs (bouts,
%  episodes, trials) for independent re-analysis in Python. Every file
%  written is listed at the end.
% =========================================================================

%% 0. PARAMETERS -----------------------------------------------------------
fs              = 200;          % frames/s
windowSec       = 115;          % analysis window after cricket drop (s).
                                % Laser runs for 120 s starting 5 s BEFORE
                                % drop, so 115 s keeps laser trials fully
                                % stimulated; applied identically to
                                % laser-off trials.
gapTolSecList   = [0 0.1 0.25]; % bridge unclassified gaps <= this inside a
                                % bout (sensitivity set)
primaryGapSec   = 0.1;          % primary analysis
tauEpSec        = 5;            % RMST horizon for episodes (s)
tauBoutSec      = 1.5;          % RMST horizon for bouts (s); with any exit as
                                % the event, bout KM curves fall to ~0 by 1.5 s
xMaxEpSec       = 10;           % x-axis limit for episode KM plots
nPerm           = 1000;         % trial-label permutations
nBoot           = 1000;         % mouse-level bootstrap resamples
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
inWinSec = accumarray(tid(inWin), 1, [nT 1]) / fs;
fprintf('  in-window frames: %.1f%% engaged, %.1f%% disengaged, %.1f%% unclassified\n', ...
    100*mean(L(inWin)==1), 100*mean(L(inWin)==-1), 100*mean(L(inWin)==0));

%% 3. BOUTS AND EPISODES FOR EACH GAP TOLERANCE; WALD HAZARD RATIOS ---------
exitNames = {'direct disengage','NOTA -> disengage','NOTA -> re-engage','capture / window end'};
lightNames = {'Light','Dark'};
R = struct();
for g = 1:numel(gapTolSecList)
    B = local_extract_bouts(L, tid, round(gapTolSecList(g)*fs), fs);
    [E, B.episode] = local_build_episodes(B);
    % drop bouts / episodes already in progress at the start of the window
    B = local_subset(B, ~B.leftTrunc);
    E = local_subset(E, ~E.leftTrunc);
    B.animal = T_animal(B.trial); B.dark = T_dark(B.trial); B.laser = T_laser(B.trial);
    E.animal = T_animal(E.trial); E.dark = T_dark(E.trial); E.laser = T_laser(E.trial);
    B.evAny = B.exitType <= 3;                    % any exit from engagement
    E.ev    = E.exitType <= 2;                    % disengagement ends the episode
    for lt = 0:1
        sE = E.dark == lt; sB = B.dark == lt;
        [b, se, p] = local_cox(E.laser(sE), E.dur(sE), E.ev(sE), E.animal(sE));
        R(g).ep.HR(lt+1) = exp(b); R(g).ep.CI(lt+1,:) = exp(b + [-1.96 1.96]*se); R(g).ep.p(lt+1) = p; R(g).ep.beta(lt+1) = b;
        [b, se, p] = local_cox(B.laser(sB), B.dur(sB), B.evAny(sB), B.animal(sB));
        R(g).bout.HR(lt+1) = exp(b); R(g).bout.CI(lt+1,:) = exp(b + [-1.96 1.96]*se); R(g).bout.p(lt+1) = p;
    end
    [~, ~, pI] = local_cox([E.laser, E.dark, E.laser.*E.dark], E.dur, E.ev, E.animal);
    R(g).ep.interactionP = pI(3);
    R(g).gapSec = gapTolSecList(g); R(g).B = B; R(g).E = E;
    fprintf('\n--- gap tolerance %.2f s: %d bouts, %d episodes ---\n', gapTolSecList(g), numel(B.dur), numel(E.dur));
    fprintf('  episodes:         laser HR light %.2f [%.2f %.2f] p=%.3g | dark %.2f [%.2f %.2f] p=%.3g | interaction p=%.3g\n', ...
        R(g).ep.HR(1), R(g).ep.CI(1,:), R(g).ep.p(1), R(g).ep.HR(2), R(g).ep.CI(2,:), R(g).ep.p(2), R(g).ep.interactionP);
    fprintf('  bouts, any exit:  laser HR light %.2f [%.2f %.2f] p=%.3g | dark %.2f [%.2f %.2f] p=%.3g\n', ...
        R(g).bout.HR(1), R(g).bout.CI(1,:), R(g).bout.p(1), R(g).bout.HR(2), R(g).bout.CI(2,:), R(g).bout.p(2));
end
gP = find(abs(gapTolSecList - primaryGapSec) < 1e-9, 1);
B = R(gP).B; E = R(gP).E;
condB = 1 + B.dark*2 + B.laser; condE = 1 + E.dark*2 + E.laser;
condNames = {'light off','light on','dark off','dark on'};

fprintf('\n=== Episodes (primary gap tolerance) ===\n');
fprintf('%12s %8s %10s %10s %10s %10s\n', '', 'n', 'median s', 'mean s', 'bouts/ep', 'P(diseng)');
for c = 1:4
    s = condE == c;
    fprintf('%12s %8d %10.2f %10.2f %10.2f %10.2f\n', condNames{c}, nnz(s), median(E.dur(s)), mean(E.dur(s)), mean(E.nBouts(s)), mean(E.ev(s)));
end

%% 4. EXIT COMPOSITION OF BOUTS (why episodes and bouts differ) -------------
exitComp = zeros(4,4);
for c = 1:4
    for e = 1:4, exitComp(c,e) = mean(B.exitType(condB==c) == e); end
end
fprintf('\n=== Bout exit composition (primary gap tolerance), rows = condition ===\n');
fprintf('%12s', ''); fprintf('%22s', exitNames{:}); fprintf('\n');
for c = 1:4, fprintf('%12s', condNames{c}); fprintf('%22.3f', exitComp(c,:)); fprintf('\n'); end

%% 5. EPISODE INFERENCE: PERMUTATION AND MOUSE BOOTSTRAP -------------------
% laser was randomized trial by trial, so trial laser labels are shuffled
% within each mouse x lighting stratum (exact null of no laser effect)
fprintf('\n=== Episode permutation test (%d perms) ===\n', nPerm);
obsB = R(gP).ep.beta;
permB = nan(nPerm, 2);
strat = T_animal*10 + T_dark;
validTrials = find(T_valid); us = unique(strat(validTrials));
for ip = 1:nPerm
    Lp = T_laser;
    for s = us'
        idx = validTrials(strat(validTrials) == s);
        Lp(idx) = T_laser(idx(randperm(numel(idx))));
    end
    lasP = Lp(E.trial);
    for lt = 0:1
        sel = E.dark == lt;
        permB(ip, lt+1) = local_cox(lasP(sel), E.dur(sel), E.ev(sel), E.animal(sel));
    end
end
pPerm = (1 + sum(abs(permB) >= abs(obsB), 1)) ./ (nPerm + 1);
fprintf('  laser HR light %.2f, p = %.3g | dark %.2f, p = %.3g\n', exp(obsB(1)), pPerm(1), exp(obsB(2)), pPerm(2));

fprintf('\n=== Episode mouse bootstrap (%d resamples) ===\n', nBoot);
mice = unique(E.animal);
bootB = nan(nBoot, 2);
for ib = 1:nBoot
    pick = mice(randi(numel(mice), numel(mice), 1));
    idx = []; newStrat = [];
    for m = 1:numel(pick)
        ii = find(E.animal == pick(m));
        idx = [idx; ii]; newStrat = [newStrat; m*ones(numel(ii),1)]; %#ok<AGROW>
    end
    for lt = 0:1
        sel = E.dark(idx) == lt;
        if numel(unique(E.laser(idx(sel)))) < 2, continue, end
        bootB(ib, lt+1) = local_cox(E.laser(idx(sel)), E.dur(idx(sel)), E.ev(idx(sel)), newStrat(sel));
    end
end
ciHR = prctile(exp(bootB), [2.5 97.5]);
ciRatio = prctile(exp(bootB(:,2) - bootB(:,1)), [2.5 97.5]);
fprintf('  laser HR light %.2f [%.2f %.2f] | dark %.2f [%.2f %.2f] | ratio dark/light %.2f [%.2f %.2f]\n', ...
    exp(obsB(1)), ciHR(:,1), exp(obsB(2)), ciHR(:,2), exp(obsB(2)-obsB(1)), ciRatio);

%% 6. PER-MOUSE RMST: EPISODES (primary) AND BOUTS (negative control) -------
rmstEp = nan(numel(mice), 4); rmstBout = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        s = E.animal == mice(m) & condE == c;
        if nnz(s) >= 5, rmstEp(m,c) = local_rmst(E.dur(s), E.ev(s), tauEpSec); end
        s = B.animal == mice(m) & condB == c;
        if nnz(s) >= 5, rmstBout(m,c) = local_rmst(B.dur(s), B.evAny(s), tauBoutSec); end
    end
end
pRm = [local_signrank(rmstEp(:,2)-rmstEp(:,1)), local_signrank(rmstEp(:,4)-rmstEp(:,3)); ...
       local_signrank(rmstBout(:,2)-rmstBout(:,1)), local_signrank(rmstBout(:,4)-rmstBout(:,3))];
fprintf('\n=== Per-mouse RMST, laser off -> on (mouse means) ===\n');
fprintf('  episodes (tau %.1f s): light %.3f -> %.3f (%d/%d up, p=%.3g) | dark %.3f -> %.3f (%d/%d up, p=%.3g)\n', tauEpSec, ...
    mean(rmstEp(:,1),'omitnan'), mean(rmstEp(:,2),'omitnan'), nnz(rmstEp(:,2)>rmstEp(:,1)), numel(mice), pRm(1,1), ...
    mean(rmstEp(:,3),'omitnan'), mean(rmstEp(:,4),'omitnan'), nnz(rmstEp(:,4)>rmstEp(:,3)), numel(mice), pRm(1,2));
fprintf('  bouts, any exit (tau %.1f s): light %.3f -> %.3f (%d/%d up, p=%.3g) | dark %.3f -> %.3f (%d/%d up, p=%.3g)\n', tauBoutSec, ...
    mean(rmstBout(:,1),'omitnan'), mean(rmstBout(:,2),'omitnan'), nnz(rmstBout(:,2)>rmstBout(:,1)), numel(mice), pRm(2,1), ...
    mean(rmstBout(:,3),'omitnan'), mean(rmstBout(:,4),'omitnan'), nnz(rmstBout(:,4)>rmstBout(:,3)), numel(mice), pRm(2,2));

%% 7. EPISODE ONSET RATE: episodes started per second outside episodes -----
nEp = accumarray(E.trial, 1, [nT 1]);
epSec = accumarray(E.trial, E.dur, [nT 1]);
outSec = inWinSec - epSec;
useT = T_valid & outSec > 1;
cT = 1 + T_dark*2 + T_laser;
onRate = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        sel = useT & T_animal == mice(m) & cT == c;
        if any(sel), onRate(m,c) = sum(nEp(sel)) / sum(outSec(sel)); end
    end
end
pOnMouse = [local_signrank(onRate(:,2)-onRate(:,1)), local_signrank(onRate(:,4)-onRate(:,3))];
RRon = nan(1,2); pOn = nan(1,2);
if exist('fitglme', 'file') && exist('table', 'file')
    TT = table(nEp(useT), log(outSec(useT)), T_laser(useT), T_dark(useT), categorical(T_animal(useT)), ...
        'VariableNames', {'n','logT','Laser','Dark','Animal'});
    for lt = 0:1
        sub = TT(TT.Dark == lt, :);
        gm = fitglme(sub, 'n ~ Laser + (1|Animal)', 'Distribution', 'Poisson', 'Offset', sub.logT);
        row = strcmp(gm.Coefficients.Name, 'Laser');
        RRon(lt+1) = exp(gm.Coefficients.Estimate(row)); pOn(lt+1) = gm.Coefficients.pValue(row);
    end
end
fprintf('\n=== Episode onset rate (per s outside episodes), laser off -> on ===\n');
fprintf('  light %.3f -> %.3f (%d/%d up, p=%.3g; GLMM RR %.2f, p=%.3g)\n', mean(onRate(:,1),'omitnan'), mean(onRate(:,2),'omitnan'), ...
    nnz(onRate(:,2)>onRate(:,1)), numel(mice), pOnMouse(1), RRon(1), pOn(1));
fprintf('  dark  %.3f -> %.3f (%d/%d up, p=%.3g; GLMM RR %.2f, p=%.3g)\n', mean(onRate(:,3),'omitnan'), mean(onRate(:,4),'omitnan'), ...
    nnz(onRate(:,4)>onRate(:,3)), numel(mice), pOnMouse(2), RRon(2), pOn(2));

%% 8. FIGURE (2 x 4) --------------------------------------------------------
if doPlot
    cOff = [0 0 0]; cOn = [0 0.6 0.9];
    fig = figure('Name', 'Engagement episodes and bouts', 'Color', 'w', 'Position', [50 100 2300 1050]);
    set(fig, 'DefaultAxesFontSize', 18);
    tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    xs = [1 2; 4 5];

    % A, B: episode survival, light and dark
    for lt = 0:1
        nexttile; hold on;
        for la = 0:1
            s = E.dark == lt & E.laser == la;
            [tt, S] = local_km(E.dur(s), E.ev(s));
            col = cOff; if la, col = cOn; end
            stairs([0; tt], [1; S], 'Color', col, 'LineWidth', 2.5);
        end
        xlim([0 xMaxEpSec]); ylim([0 1]); box off;
        xlabel('Episode duration (s)'); ylabel('P(still engaged)');
        title(sprintf('%c. %s episodes: HR %.2f, p = %.2g', 'A'+lt, lightNames{lt+1}, exp(obsB(lt+1)), pPerm(lt+1)));
        legend({sprintf('laser off (%d)', nnz(E.dark==lt & E.laser==0)), sprintf('laser on (%d)', nnz(E.dark==lt & E.laser==1))}, ...
            'Location', 'northeast', 'Box', 'off');
    end

    % C: continuous bouts (any exit), light dashed, dark solid
    nexttile; hold on;
    sty = {'--','-'};
    for lt = 0:1
        for la = 0:1
            s = B.dark == lt & B.laser == la;
            [tt, S] = local_km(B.dur(s), B.evAny(s));
            col = cOff; if la, col = cOn; end
            stairs([0; tt], [1; S], sty{lt+1}, 'Color', col, 'LineWidth', 2.5);
        end
    end
    xlim([0 tauBoutSec]); ylim([0 1]); box off;
    xlabel('Bout duration (s)'); ylabel('P(still in bout)');
    legend(condNames, 'Location', 'northeast', 'Box', 'off');
    title(sprintf('C. Continuous bouts: HR light %.2f, dark %.2f', R(gP).bout.HR));

    % D: sensitivity forest plot
    nexttile; hold on;
    y = 0; yv = []; yl = {};
    for lt = 0:1
        col = cOff; if lt, col = cOn; end
        for kind = 1:2
            for g = 1:numel(gapTolSecList)
                y = y + 1;
                if kind == 1, hr = R(g).ep.HR(lt+1); ci = R(g).ep.CI(lt+1,:); mk = 'o'; lab = 'episode';
                else,         hr = R(g).bout.HR(lt+1); ci = R(g).bout.CI(lt+1,:); mk = 's'; lab = 'bout'; end
                plot(ci, [y y], '-', 'Color', col, 'LineWidth', 2);
                plot(hr, y, mk, 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 9);
                tag = ''; if g == gP && kind == 1, tag = ' *'; end
                yv(end+1) = y; yl{end+1} = sprintf('%s %s, gap %.2f%s', lightNames{lt+1}, lab, gapTolSecList(g), tag); %#ok<SAGROW>
            end
        end
        y = y + 1;
    end
    xline(1, '--k'); set(gca, 'XScale', 'log', 'YTick', yv, 'YTickLabel', yl, 'YDir', 'reverse');
    xlabel('Laser hazard ratio (< 1 = longer)'); title('D. Sensitivity (Wald 95% CI)'); box off;

    % E, F: per-mouse RMST (episodes, bouts)
    rm = {rmstEp, rmstBout}; taus = [tauEpSec tauBoutSec]; nm = {'episodes','bouts'};
    for k = 1:2
        nexttile; hold on;
        for lt = 0:1
            cc = lt*2 + (1:2);
            plot(xs(lt+1,:), rm{k}(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
            plot(xs(lt+1,1), rm{k}(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
            plot(xs(lt+1,2), rm{k}(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
        end
        set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]); box off;
        local_lightdark_labels();
        ylabel(sprintf('RMST to %.1f s (s)', taus(k)));
        title(sprintf('%c. Per mouse, %s: p = %.2g, %.2g', 'D'+k, nm{k}, pRm(k,1), pRm(k,2)));
    end

    % G: episode onset rate
    nexttile; hold on;
    for lt = 0:1
        cc = lt*2 + (1:2);
        plot(xs(lt+1,:), onRate(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(xs(lt+1,1), onRate(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
        plot(xs(lt+1,2), onRate(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
    end
    set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]); box off;
    ylabel('Episodes / s outside episodes');
    local_lightdark_labels();
    title(sprintf('G. Episode onsets: p = %.2g, %.2g', pOnMouse)); box off;

    % H: bout exit composition
    nexttile;
    bh = bar(exitComp, 'stacked');
    exitCols = [0.80 0.20 0.20; 0.95 0.60 0.50; 0.60 0.60 0.60; 0.30 0.60 0.30];
    for e = 1:4, bh(e).FaceColor = exitCols(e,:); end
    set(gca, 'XTickLabel', condNames, 'XTickLabelRotation', 30); ylabel('Fraction of bouts'); ylim([0 1]); box off;
    legend(exitNames, 'Location', 'southoutside', 'Box', 'off');
    title('H. How bouts end');

    title(tl, sprintf(['Engagement episodes (primary) and continuous bouts; gap tolerance %.2f s, first %d s after drop. ' ...
        'p-values: A, B permutation; E-G light, dark signed-rank'], primaryGapSec, windowSec), 'FontSize', 18);
    outFile = fullfile(outputrootdir, 'engaged_bout_survival.pdf'); exportgraphics(fig, outFile, 'ContentType', 'vector'); written{end+1} = outFile;
end

%% 9. SAVE RESULTS AND EXPORT TABLES FOR INDEPENDENT CHECKS -----------------
Rsave = rmfield(R, {'B','E'});
results = struct('R', Rsave, 'pPerm', pPerm, 'ciHR', ciHR, 'ciRatio', ciRatio, 'rmstEp', rmstEp, ...
    'rmstBout', rmstBout, 'pRm', pRm, 'onRate', onRate, 'pOnMouse', pOnMouse, 'RRon', RRon, 'pOn', pOn, ...
    'exitComp', exitComp, 'mice', mice, 'params', struct('fs', fs, 'windowSec', windowSec, ...
    'gapTolSecList', gapTolSecList, 'primaryGapSec', primaryGapSec, 'tauEpSec', tauEpSec, ...
    'tauBoutSec', tauBoutSec, 'nPerm', nPerm, 'nBoot', nBoot));
outFile = fullfile(outputrootdir, 'engaged_bout_survival_results.mat'); save(outFile, 'results'); written{end+1} = outFile;
if exist('writetable', 'file') && exist('table', 'file')
    BT = table(B.trial, uf(B.trial), B.animal, B.dark, B.laser, B.startFrame, B.dur, B.exitType, B.episode, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','start_frame','dur_s','exit_type','episode_id'});
    outFile = fullfile(outputrootdir, 'engaged_bouts_primary.csv'); writetable(BT, outFile); written{end+1} = outFile;
    ET = table(E.id, E.trial, uf(E.trial), E.animal, E.dark, E.laser, E.startFrame, E.dur, E.nBouts, E.exitType, E.ev, ...
        'VariableNames', {'episode_id','trial_id','filename','animal','dark','laser','start_frame','dur_s','n_bouts','last_exit_type','event'});
    outFile = fullfile(outputrootdir, 'engagement_episodes_primary.csv'); writetable(ET, outFile); written{end+1} = outFile;
    TTr = table((1:nT)', uf, T_animal, T_dark, T_laser, T_drop, inWinSec, nEp, epSec, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','cricketdrop','inwindow_s','n_episodes','episode_s'});
    outFile = fullfile(outputrootdir, 'engaged_bout_trials_primary.csv'); writetable(TTr, outFile); written{end+1} = outFile;
end
fprintf('\nFiles written:\n'); fprintf('  %s\n', written{:}); fprintf('Done.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function B = local_extract_bouts(L, tid, gapFrames, fs)
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
B.endFrame = re(ib);
B.dur = (re(ib) - rs(ib) + 1) / fs;
B.exitType = exitType;
B.leftTrunc = leftTrunc;
B.fs = fs;
end

function [E, epId] = local_build_episodes(B)
% Chain consecutive bouts within a trial across exit type 3 (NOTA -> re-engage).
newEp = [true; B.exitType(1:end-1) ~= 3 | diff(B.trial) ~= 0];
epId = cumsum(newEp);
nE = epId(end);
firstIdx = find(newEp);
lastIdx = [firstIdx(2:end) - 1; numel(epId)];
E.id = (1:nE)';
E.trial = B.trial(firstIdx);
E.startFrame = B.startFrame(firstIdx);
E.dur = (B.endFrame(lastIdx) - B.startFrame(firstIdx) + 1) / B.fs;
E.nBouts = lastIdx - firstIdx + 1;
E.exitType = B.exitType(lastIdx);
E.leftTrunc = B.leftTrunc(firstIdx);
end

function S = local_subset(S, keep)
fn = fieldnames(S);
for f = 1:numel(fn)
    if numel(S.(fn{f})) == numel(keep), S.(fn{f}) = S.(fn{f})(keep); end
end
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

function p = local_signrank(d)
d = d(~isnan(d));
if numel(d) < 2, p = NaN; return, end
p = signrank(d);
end

function local_lightdark_labels()
% 'light' / 'dark' group labels inside the top of a paired off/on panel
yl = ylim; ylim(yl);
text(1.5, yl(2), 'light', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 18);
text(4.5, yl(2), 'dark',  'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 18);
end
