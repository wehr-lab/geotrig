%% =========================================================================
%  ANALYZE_REENGAGE_GAPS
%  What happens in the gap after an engaged bout ends? Tests whether
%  "NOTA -> re-engage" exits are contacts (or near-contacts), and whether
%  ZIm activation increases re-pursuit after a contact.
%
%  Companion to analyze_engaged_bout_survival.m: same frame labels, same
%  analysis window, same gap tolerance, same exit types.
%
%  For every engaged bout that does not run into capture / window end:
%    exit 1 (direct disengage)  gap = first 0.25 s of the disengaged state
%    exit 2 (NOTA -> disengage) gap = the NOTA run
%    exit 3 (NOTA -> re-engage) gap = the NOTA run
%  Geometry is computed on 0.25 s median-filtered signals (the window the
%  follow/chase/hot pursuit detectors use).
%
%  Gap category (first rule that applies, in this order):
%    1 contact      min range < 5 cm  (intercept threshold; with hot pursuit
%                   at 5-20 cm, every frame < 5 cm is outside all engaged states)
%    2 far approach median range > 20 cm and median |az| < 30 deg
%    3 fast off-axis median range 5-20 cm, |az| > 30 deg, speed > 15 cm/s
%    4 still cricket stalk geometry but median cricket speed < 1 cm/s
%    5 other
%
%  INPUTS (workspace, Nframes x 1): everything analyze_engaged_bout_survival.m
%  needs, plus range, az, mouse_spd, cricket_spd. Optional (absolute frame
%  indices, used if present): cricket_jump_event_frames,
%  contact_gain_event_frames, intercept_event_frames.
%
%  OUTPUTS: figure (6 panels), results struct, gap table CSV (one row per
%  gap) and trial table CSV, for independent re-analysis.
% =========================================================================

%% 0. PARAMETERS -----------------------------------------------------------
fs            = 200;
windowSec     = 115;      % same as analyze_engaged_bout_survival.m
gapTolSec     = 0.1;      % primary gap tolerance
medWin        = round(0.25*fs);
contactCm     = 5;        % contact = min range below this
nearCm        = 10;       % detect_contact.m threshold, for the pause check
evPreSec      = 0.25;     % events this long before a gap also count for it
snipPre = 1; snipPost = 2;% seconds around gap onset for aligned traces
nPerm         = 1000;
doPlot        = true;
rng(1);

need = {'hotpursuit','chase','follow','stalk','wander','pause','cricket_present', ...
        'localframe','filename','metadata','range','az','mouse_spd','cricket_spd'};
hasVar = false(size(need));
for iv = 1:numel(need), hasVar(iv) = exist(need{iv}, 'var') == 1; end
missing = need(~hasVar);
if ~isempty(missing), error('Missing inputs: %s', strjoin(missing, ', ')); end
if ~exist('outputrootdir','var') || ~exist(outputrootdir,'dir'), outputrootdir = pwd; end
written = {};                   % every file this script writes, listed at the end

%% 1. TRIAL INFO (identical to analyze_engaged_bout_survival.m) ------------
fnames = cellstr(filename(:));
[uf, ~, tid] = unique(fnames, 'stable');
tid = tid(:); nT = numel(uf);
mdFile = cellstr(metadata.filename);
T_drop = nan(nT,1); T_laser = nan(nT,1); T_dark = nan(nT,1); T_animal = nan(nT,1);
for k = 1:nT
    r = find(cellfun(@(m) ~isempty(strfind(uf{k}, m)), mdFile), 1);
    if isempty(r), continue, end
    T_drop(k) = metadata.cricketdrop(r);
    lv = metadata.laser_value(r); if isnan(lv), lv = 0; end
    T_laser(k) = double(lv ~= 0);
    T_dark(k) = double(metadata.condition(r) == 2);
    tok = regexp(uf{k}, 'mouse-(\d+)', 'tokens', 'once');
    if ~isempty(tok), T_animal(k) = str2double(tok{1}); end
end
T_valid = ~isnan(T_drop) & T_drop > 0 & ~isnan(T_animal) & ~isnan(T_dark);

%% 2. FRAME LABELS AND STATE CODES -----------------------------------------
eng = logical(hotpursuit(:)) | logical(chase(:)) | logical(follow(:));
dis = (logical(stalk(:)) | logical(wander(:)) | logical(pause(:))) & ~eng;
tsd = (localframe(:) - T_drop(tid)) / fs;
inWin = logical(cricket_present(:)) & tsd >= 0 & tsd <= windowSec & T_valid(tid);
L = nan(numel(eng),1); L(inWin) = 0; L(inWin & dis) = -1; L(inWin & eng) = 1;
S = zeros(numel(eng),1);                         % 1 hp 2 ch 3 fo 4 st 5 wa 6 pa
S(logical(pause(:)) & dis) = 6; S(logical(wander(:)) & dis) = 5; S(logical(stalk(:)) & dis) = 4;
S(logical(follow(:))) = 3; S(logical(chase(:))) = 2; S(logical(hotpursuit(:))) = 1;
stateNames = {'hot pursuit','chase','follow','stalk','wander','pause'};

mR = medfilt1(range(:),       medWin, 'omitnan');
mA = medfilt1(abs(az(:)),     medWin, 'omitnan');
mV = medfilt1(mouse_spd(:),   medWin, 'omitnan');
mC = medfilt1(cricket_spd(:), medWin, 'omitnan');

%% 3. EXTRACT GAPS ---------------------------------------------------------
G = local_extract_gaps(L, S, tid, round(gapTolSec*fs), fs);
nG = numel(G.trial);
fprintf('\n=== %d gaps after engaged bouts (exit 1: %d, exit 2: %d, exit 3: %d) ===\n', ...
    nG, nnz(G.exitType==1), nnz(G.exitType==2), nnz(G.exitType==3));

G.animal = T_animal(G.trial); G.dark = T_dark(G.trial); G.laser = T_laser(G.trial);
G.dur = (G.gapEnd - G.gapStart + 1) / fs;
[G.minRange, G.medRange, G.frac5, G.frac10, G.medSpd, G.medAz, G.medCr, G.maxCr] = deal(nan(nG,1));
for k = 1:nG
    idx = G.gapStart(k):G.gapEnd(k);
    G.minRange(k) = min(mR(idx));           G.medRange(k) = median(mR(idx), 'omitnan');
    G.frac5(k)    = mean(mR(idx) < contactCm); G.frac10(k) = mean(mR(idx) < nearCm);
    G.medSpd(k)   = median(mV(idx), 'omitnan'); G.medAz(k)  = median(mA(idx), 'omitnan');
    G.medCr(k)    = median(mC(idx), 'omitnan'); G.maxCr(k)  = max(cricket_spd(idx));
end
G.rangeAtNext = mR(G.nextStart);

% category (assigned lowest priority first, so higher-priority rules win)
G.category = 5*ones(nG,1);
G.category(G.medSpd > 2 & G.medSpd < 15 & G.medRange > 5 & G.medRange < 20 & G.medAz > 30 & G.medCr < 1) = 4;
G.category(G.medAz > 30 & G.medSpd > 15 & G.medRange >= 5 & G.medRange <= 20) = 3;
G.category(G.medRange > 20 & G.medAz < 30) = 2;
G.category(G.minRange < contactCm) = 1;
catNames = {'contact (<5 cm)','far approach','fast off-axis','still cricket','other'};

% handling check: disengagements that land in pause right next to the cricket
G.pauseNear = ismember(G.exitType, [1 2]) & G.nextState == 6 & G.rangeAtNext < nearCm;

% optional event overlaps
evList = {'cricket_jump_event_frames','contact_gain_event_frames','intercept_event_frames'};
evShort = {'evJump','evContact','evIntercept'};
for e = 1:numel(evList)
    G.(evShort{e}) = nan(nG,1);
    if exist(evList{e}, 'var')
        ev = sort(eval(evList{e})); ev = ev(:);
        lo = G.gapStart - round(evPreSec*fs); hi = G.gapEnd;
        G.(evShort{e}) = double(arrayfun(@(a,b) any(ev >= a & ev <= b), lo, hi));
    end
end

% aligned range traces around gap onset (NaN outside the trial's window)
tSnip = (-round(snipPre*fs):round(snipPost*fs)) / fs;
snipIdx = G.gapStart + round(tSnip*fs);
snipIdx(snipIdx < 1 | snipIdx > numel(mR)) = 1;
okIdx = inWin(snipIdx) & tid(snipIdx) == G.trial;
snipR = mR(snipIdx); snipR(~okIdx) = NaN;

%% 4. Q1: ARE RE-ENGAGE GAPS CONTACTS? ------------------------------------
condIdx = 1 + G.dark*2 + G.laser;
condNames = {'light off','light on','dark off','dark on'};
outNames = {'direct disengage','NOTA -> disengage','NOTA -> re-engage'};
fprintf('\n=== Q1. Gap category by exit type (all conditions) ===\n');
fprintf('%20s', ''); fprintf('%17s', catNames{:}); fprintf('%8s\n', 'n');
comp = zeros(3,5);
for et = 1:3
    for c = 1:5, comp(et,c) = mean(G.category(G.exitType==et) == c); end
    fprintf('%20s', outNames{et}); fprintf('%17.3f', comp(et,:)); fprintf('%8d\n', nnz(G.exitType==et));
end
fprintf('  median min range: re-engage %.1f cm, NOTA->disengage %.1f cm, direct %.1f cm\n', ...
    median(G.minRange(G.exitType==3),'omitnan'), median(G.minRange(G.exitType==2),'omitnan'), median(G.minRange(G.exitType==1),'omitnan'));
fprintf('  disengagements into pause within %d cm of the cricket: %d of %d (%.1f%%)\n', ...
    nearCm, nnz(G.pauseNear), nnz(G.exitType<=2), 100*nnz(G.pauseNear)/nnz(G.exitType<=2));
compCond = zeros(4,5,2);                 % condition x category x {re-engage, disengage}
for c = 1:4
    for k = 1:5
        compCond(c,k,1) = mean(G.category(condIdx==c & G.exitType==3) == k);
        compCond(c,k,2) = mean(G.category(condIdx==c & G.exitType<=2) == k);
    end
end

%% 5. Q2: DOES THE LASER INCREASE RE-PURSUIT AFTER A CONTACT? --------------
isC = G.category == 1;
reeng = G.exitType == 3;
mice = unique(G.animal(~isnan(G.animal)));
pRe = nan(numel(mice), 4); nRe = zeros(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        sel = isC & G.animal == mice(m) & condIdx == c;
        nRe(m,c) = nnz(sel);
        if nnz(sel) >= 5, pRe(m,c) = mean(reeng(sel)); end
    end
end
obsOR = [local_mh_logor(reeng(isC & G.dark==0), G.laser(isC & G.dark==0), G.animal(isC & G.dark==0)), ...
         local_mh_logor(reeng(isC & G.dark==1), G.laser(isC & G.dark==1), G.animal(isC & G.dark==1))];
% permutation: trial laser labels shuffled within mouse x lighting
permOR = nan(nPerm, 2);
validTrials = find(T_valid); strat = T_animal*10 + T_dark; us = unique(strat(validTrials));
for ip = 1:nPerm
    Lp = T_laser;
    for s = us'
        idx = validTrials(strat(validTrials) == s);
        Lp(idx) = T_laser(idx(randperm(numel(idx))));
    end
    lp = Lp(G.trial);
    for lt = 0:1
        sel = isC & G.dark == lt;
        permOR(ip, lt+1) = local_mh_logor(reeng(sel), lp(sel), G.animal(sel));
    end
end
pPermOR = (1 + sum(abs(permOR) >= abs(obsOR), 1)) ./ (nPerm + 1);
[~, pMouseL] = local_signrank(pRe(:,2) - pRe(:,1));
[~, pMouseD] = local_signrank(pRe(:,4) - pRe(:,3));
fprintf('\n=== Q2. P(re-engage | contact), laser off -> on ===\n');
fprintf('  light: %.2f -> %.2f (mouse means), MH odds ratio %.2f, permutation p = %.3g, %d/%d mice up (signed-rank p = %.3g)\n', ...
    mean(pRe(:,1),'omitnan'), mean(pRe(:,2),'omitnan'), exp(obsOR(1)), pPermOR(1), nnz(pRe(:,2) > pRe(:,1)), nnz(all(~isnan(pRe(:,1:2)),2)), pMouseL);
fprintf('  dark:  %.2f -> %.2f (mouse means), MH odds ratio %.2f, permutation p = %.3g, %d/%d mice up (signed-rank p = %.3g)\n', ...
    mean(pRe(:,3),'omitnan'), mean(pRe(:,4),'omitnan'), exp(obsOR(2)), pPermOR(2), nnz(pRe(:,4) > pRe(:,3)), nnz(all(~isnan(pRe(:,3:4)),2)), pMouseD);

glmeOut = struct();
if exist('fitglme','file') && exist('table','file')
    crC = G.maxCr(isC);
    crZ = (crC - mean(crC, 'omitnan')) ./ std(crC, 'omitnan');   % NaN-safe z-score (zscore returns all NaN if any are NaN)
    fprintf('  contact gaps with no cricket-speed samples (dropped from the covariate model only): %d\n', nnz(isnan(crZ)));
    GT = table(double(reeng(isC)), G.laser(isC), G.dark(isC), categorical(G.animal(isC)), ...
        categorical(G.trial(isC)), crZ, 'VariableNames', {'Reeng','Laser','Dark','Animal','Trial','CrMax'});
    gmI = fitglme(GT, 'Reeng ~ Laser*Dark + (1|Animal) + (1|Trial)', 'Distribution', 'Binomial');
    glmeOut.interaction = gmI.Coefficients(strcmp(gmI.Coefficients.Name, 'Laser:Dark'), :);
    for lt = 0:1
        sub = GT(GT.Dark == lt, :);
        g1 = fitglme(sub, 'Reeng ~ Laser + (1|Animal) + (1|Trial)', 'Distribution', 'Binomial');
        r1 = strcmp(g1.Coefficients.Name, 'Laser');
        glmeOut.OR(lt+1) = exp(g1.Coefficients.Estimate(r1)); glmeOut.p(lt+1) = g1.Coefficients.pValue(r1);
        nCr = nnz(~isnan(sub.CrMax));
        ltName = {'light','dark'};
        fprintf('  %s: %d of %d contact gaps have cricket-speed samples\n', ltName{lt+1}, nCr, height(sub));
        glmeOut.ORcr(lt+1) = NaN; glmeOut.pcr(lt+1) = NaN;
        if nCr >= 50 && numel(unique(sub.Laser(~isnan(sub.CrMax)))) == 2
            g2 = fitglme(sub(~isnan(sub.CrMax),:), 'Reeng ~ Laser + CrMax + (1|Animal) + (1|Trial)', 'Distribution', 'Binomial');
            r2 = strcmp(g2.Coefficients.Name, 'Laser');
            glmeOut.ORcr(lt+1) = exp(g2.Coefficients.Estimate(r2)); glmeOut.pcr(lt+1) = g2.Coefficients.pValue(r2);
        else
            fprintf('    too few contact gaps with cricket speed; covariate model skipped\n');
        end
    end
    fprintf('  binomial GLMM, laser OR: light %.2f (p = %.3g), dark %.2f (p = %.3g); laser x dark p = %.3g\n', ...
        glmeOut.OR(1), glmeOut.p(1), glmeOut.OR(2), glmeOut.p(2), glmeOut.interaction.pValue);
    fprintf('  controlling for peak cricket speed: light %.2f (p = %.3g), dark %.2f (p = %.3g)\n', ...
        glmeOut.ORcr(1), glmeOut.pcr(1), glmeOut.ORcr(2), glmeOut.pcr(2));
end

%% 6. Q3: DOES THE LASER CHANGE HOW OFTEN CONTACTS HAPPEN? ------------------
engSec = accumarray(tid(L == 1), 1, [nT 1]) / fs;
nCont  = accumarray(G.trial(isC), 1, [nT 1]);
cT = 1 + T_dark*2 + T_laser;
cRate = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        sel = T_valid & T_animal == mice(m) & cT == c & engSec > 0;
        if any(sel), cRate(m,c) = sum(nCont(sel)) / sum(engSec(sel)); end
    end
end
[~, pCrL] = local_signrank(cRate(:,2) - cRate(:,1));
[~, pCrD] = local_signrank(cRate(:,4) - cRate(:,3));
fprintf('\n=== Q3. Contacts per second engaged, laser off -> on ===\n');
fprintf('  light: %.3f -> %.3f (signed-rank p = %.3g)\n', mean(cRate(:,1),'omitnan'), mean(cRate(:,2),'omitnan'), pCrL);
fprintf('  dark:  %.3f -> %.3f (signed-rank p = %.3g)\n', mean(cRate(:,3),'omitnan'), mean(cRate(:,4),'omitnan'), pCrD);

%% 7. Q4: HOW FAST DO MICE RE-ENGAGE AFTER A CONTACT? -----------------------
% cumulative fraction of contact gaps re-engaged by time t (disengagement
% competes; it is not censoring, so this is a cumulative incidence)
tGrid = 0:0.05:3;
cumRe = nan(4, numel(tGrid));
for c = 1:4
    sel = isC & condIdx == c;
    cumRe(c,:) = mean(reeng(sel) & G.dur(sel) <= tGrid, 1);
end
medLatRe = nan(numel(mice), 4);
for m = 1:numel(mice)
    for c = 1:4
        sel = isC & reeng & G.animal == mice(m) & condIdx == c;
        if nnz(sel) >= 5, medLatRe(m,c) = median(G.dur(sel)); end
    end
end
fprintf('\n=== Q4. Median contact -> re-engage latency (s), mouse means, off -> on ===\n');
fprintf('  light: %.2f -> %.2f   dark: %.2f -> %.2f\n', mean(medLatRe,1,'omitnan'));

%% 8. Q5: CRICKET ESCAPE CONTROL -------------------------------------------
fprintf('\n=== Q5. Peak cricket speed during contact gaps (cm/s), median off -> on ===\n');
for lt = 0:1
    a = G.maxCr(isC & G.dark==lt & G.laser==0); b = G.maxCr(isC & G.dark==lt & G.laser==1);
    fprintf('  %s: %.1f -> %.1f', condNames{lt*2+1}(1:5), median(a,'omitnan'), median(b,'omitnan'));
    if ~all(isnan(G.evJump)), fprintf(';  P(cricket jump) %.2f -> %.2f', mean(G.evJump(isC & G.dark==lt & G.laser==0)), mean(G.evJump(isC & G.dark==lt & G.laser==1))); end
    fprintf('\n');
end

%% 9. FIGURE ---------------------------------------------------------------
if doPlot
    cOff = [0 0 0]; cOn = [0 0.6 0.9];
    catCols = [0.85 0.33 0.10; 0.30 0.60 0.30; 0.55 0.40 0.75; 0.95 0.75 0.20; 0.65 0.65 0.65];
    fig = figure('Name','Re-engagement gaps','Color','w','Position',[100 100 1800 1050]);
    set(fig, 'DefaultAxesFontSize', 18);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % A: minimum range during the gap, by exit type
    nexttile; hold on;
    edges = 0:1:40; ctr = edges(1:end-1) + 0.5;
    exCols = [0.80 0.20 0.20; 0.95 0.60 0.50; 0.45 0.45 0.45];
    for et = 1:3
        h = histcounts(G.minRange(G.exitType==et), edges, 'Normalization', 'probability');
        plot(ctr, h, '-', 'Color', exCols(et,:), 'LineWidth', 2.5);
    end
    xline(contactCm, '--k'); xline(nearCm, ':k');
    xlabel('Min range during gap (cm)'); ylabel('Fraction of gaps');
    legend(outNames, 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
    title('A. Are re-engage gaps contacts?'); box off;

    % B: category composition, re-engage vs disengage gaps, per condition
    nexttile;
    M = [compCond(:,:,1); compCond(:,:,2)];
    bh = bar(M, 'stacked');
    for k = 1:5, bh(k).FaceColor = catCols(k,:); end
    set(gca, 'XTick', 1:8, 'XTickLabel', [strcat('R:', condNames), strcat('D:', condNames)], 'XTickLabelRotation', 45);
    ylabel('Fraction of gaps'); ylim([0 1]);
    legend(catNames, 'Location', 'southoutside', 'Box', 'off', 'FontSize', 18, 'NumColumns', 2);
    title('B. Gap type (R = re-engage, D = disengage)'); box off;

    % C: per-mouse P(re-engage | contact)
    nexttile; hold on; xs = [1 2; 4 5];
    for lt = 0:1
        cc = lt*2 + (1:2);
        plot(xs(lt+1,:), pRe(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(xs(lt+1,1), pRe(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
        plot(xs(lt+1,2), pRe(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
    end
    set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]); ylim([0 1]);
    text(1.5, 0.05, 'light', 'HorizontalAlignment', 'center', 'FontSize', 18);
    text(4.5, 0.05, 'dark',  'HorizontalAlignment', 'center', 'FontSize', 18);
    ylabel('P(re-engage | contact)');
    title(sprintf('C. OR light %.2f (p=%.2g), dark %.2f (p=%.2g)', exp(obsOR(1)), pPermOR(1), exp(obsOR(2)), pPermOR(2))); box off;

    % D: per-mouse contact rate
    nexttile; hold on;
    for lt = 0:1
        cc = lt*2 + (1:2);
        plot(xs(lt+1,:), cRate(:,cc)', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(xs(lt+1,1), cRate(:,cc(1)), 'o', 'MarkerFaceColor', cOff, 'MarkerEdgeColor', cOff, 'MarkerSize', 9);
        plot(xs(lt+1,2), cRate(:,cc(2)), 'o', 'MarkerFaceColor', cOn,  'MarkerEdgeColor', cOn,  'MarkerSize', 9);
    end
    set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'off','on','off','on'}); xlim([0.3 5.7]);
    ylabel('Contacts / s engaged');
    title(sprintf('D. Contact rate: light p=%.2g, dark p=%.2g', pCrL, pCrD)); box off;

    % E: cumulative re-engagement after contact
    nexttile; hold on;
    sty = {'--','--','-','-'}; cols = {cOff, cOn, cOff, cOn};
    for c = 1:4, plot(tGrid, cumRe(c,:), sty{c}, 'Color', cols{c}, 'LineWidth', 2.5); end
    xlabel('Time since contact gap onset (s)'); ylabel('Fraction re-engaged');
    legend(condNames, 'Location', 'southeast', 'Box', 'off', 'FontSize', 18);
    title('E. Re-engagement after contact'); box off;

    % F: aligned range around contact gaps, dark, by outcome and laser
    nexttile; hold on;
    for la = 0:1
        for oc = [1 0]
            sel = isC & G.dark == 1 & G.laser == la & reeng == oc;
            col = cOff; if la, col = cOn; end
            ls = '-'; if ~oc, ls = ':'; end
            plot(tSnip, mean(snipR(sel,:), 1, 'omitnan'), ls, 'Color', col, 'LineWidth', 2.5);
        end
    end
    xline(0, '--k'); yline(contactCm, ':k');
    xlabel('Time from gap onset (s)'); ylabel('Range (cm)');
    legend({'off, re-engage','off, disengage','on, re-engage','on, disengage'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
    title('F. Dark: range around contact gaps'); box off;

    title(tl, sprintf('Gaps after engaged bouts (gap tolerance %.2f s, first %d s after drop)', gapTolSec, windowSec), 'FontSize', 18);
    outFile = fullfile(outputrootdir, 'reengage_gaps.pdf'); exportgraphics(fig, outFile, 'ContentType', 'vector'); written{end+1} = outFile;
end

%% 10. SAVE AND EXPORT ------------------------------------------------------
results = struct('comp', comp, 'compCond', compCond, 'pRe', pRe, 'nRe', nRe, 'obsOR', exp(obsOR), ...
    'pPermOR', pPermOR, 'glme', glmeOut, 'cRate', cRate, 'medLatRe', medLatRe, 'cumRe', cumRe, ...
    'tGrid', tGrid, 'mice', mice, 'params', struct('fs', fs, 'windowSec', windowSec, 'gapTolSec', gapTolSec, ...
    'medWin', medWin, 'contactCm', contactCm, 'nearCm', nearCm, 'nPerm', nPerm));
outFile = fullfile(outputrootdir, 'reengage_gaps_results.mat'); save(outFile, 'results', 'tSnip', 'snipR'); written{end+1} = outFile;
if exist('writetable','file') && exist('table','file')
    GTab = table(G.trial, uf(G.trial), G.animal, G.dark, G.laser, G.exitType, G.nextState, G.gapStart, G.gapEnd, ...
        G.dur, G.minRange, G.medRange, G.frac5, G.frac10, G.medSpd, G.medAz, G.medCr, G.maxCr, G.category, ...
        G.pauseNear, G.evJump, G.evContact, G.evIntercept, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','exit_type','next_state','gap_start','gap_end', ...
        'dur_s','min_range','med_range','frac_lt5','frac_lt10','med_spd','med_az','med_cricket_spd','max_cricket_spd', ...
        'category','pause_near','ev_jump','ev_contact','ev_intercept'});
    outFile = fullfile(outputrootdir, 'reengage_gaps.csv'); writetable(GTab, outFile); written{end+1} = outFile;
    TTr = table((1:nT)', uf, T_animal, T_dark, T_laser, engSec, nCont, ...
        'VariableNames', {'trial_id','filename','animal','dark','laser','engaged_s','n_contact_gaps'});
    outFile = fullfile(outputrootdir, 'reengage_gap_trials.csv'); writetable(TTr, outFile); written{end+1} = outFile;
end
fprintf('\nFiles written:\n'); fprintf('  %s\n', written{:}); fprintf('Done.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function G = local_extract_gaps(L, S, tid, gapFrames, fs)
% Same run-length encoding and gap bridging as analyze_engaged_bout_survival.m;
% returns the gap after each engaged bout that ends in exit type 1, 2 or 3.
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
nr = numel(rs); ib = find(rv == 1);
out = nan(numel(ib), 7); k = 0;
for j = 1:numel(ib)
    i = ib(j);
    if i == 1 || rt(i-1) ~= rt(i) || rv(i-1) == -9, continue, end   % engaged at window start
    if i == nr || rt(i+1) ~= rt(i) || rv(i+1) == -9, continue, end  % capture / window end
    if rv(i+1) == -1
        et = 1; gs = rs(i+1); ge = min(re(i+1), rs(i+1) + round(0.25*fs) - 1); ns = rs(i+1);
    elseif i+2 <= nr && rt(i+2) == rt(i) && rv(i+2) == -1
        et = 2; gs = rs(i+1); ge = re(i+1); ns = rs(i+2);
    elseif i+2 <= nr && rt(i+2) == rt(i) && rv(i+2) == 1
        et = 3; gs = rs(i+1); ge = re(i+1); ns = rs(i+2);
    else
        continue                                                     % NOTA then capture / window end
    end
    k = k + 1;
    out(k,:) = [rt(i), rs(i), re(i), gs, ge, et, ns];
end
out = out(1:k,:);
G.trial = out(:,1); G.boutStart = out(:,2); G.boutEnd = out(:,3);
G.gapStart = out(:,4); G.gapEnd = out(:,5); G.exitType = out(:,6);
G.nextStart = out(:,7); G.nextState = S(out(:,7));
end

function [rs, re, rv, rt] = local_runs(v, tid)
n = numel(v);
chg = [true; diff(v) ~= 0 | diff(tid) ~= 0];
rs = find(chg); re = [rs(2:end) - 1; n];
rv = v(rs); rt = tid(rs);
end

function lor = local_mh_logor(y, x, strata)
% Mantel-Haenszel common log odds ratio of y (0/1) on x (0/1) across strata
y = logical(y); x = logical(x);
us = unique(strata); num = 0; den = 0;
for s = us'
    k = strata == s; n = nnz(k);
    if n == 0, continue, end
    a = nnz(k & x & y); b = nnz(k & x & ~y); c = nnz(k & ~x & y); d = nnz(k & ~x & ~y);
    num = num + a*d/n; den = den + b*c/n;
end
lor = log((num + 0.5) / (den + 0.5));
end

function [h, p] = local_signrank(d)
d = d(~isnan(d));
if numel(d) < 2, h = false; p = NaN; return, end
p = signrank(d); h = p < 0.05;
end
