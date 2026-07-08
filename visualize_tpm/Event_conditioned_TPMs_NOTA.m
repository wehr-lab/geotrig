% Event_conditioned_TPMs_NOTA

% Event_conditioned_TPMs  (revised)
%
% Changes from original:
%   1. NOTA (none-of-the-above) state: frames where no labelled state is
%      active are assigned state index (nStates+1) rather than 0.  They are
%      carried through the sequence but transitions that touch NOTA are
%      excluded from the *partial* TPM counts — so NOTA acts as a censoring
%      label rather than a true state.  This keeps the frame timeline intact
%      (important for event-window bookkeeping) while restricting inference
%      to labelled-state transitions only.
%
%   2. Partial TPM: only transitions where BOTH from- and to-state are
%      labelled (1..nStates) are counted.  The resulting nStates x nStates
%      TPM is row-normalised over labelled destinations only.
%
%   3. Condition masking (optional, commented block): applying a condition
%      mask BEFORE building stateSeq is the cleanest approach — just AND
%      each column of stateMask with the condition vector.  Frames outside
%      the condition become NOTA and are automatically censored.  This means
%      the event-window logic is unaffected; windows that straddle a
%      condition boundary will simply have fewer valid transitions inside
%      them, which is the correct behaviour.


% What changed and why 
% NOTA state (NOTA = nStates+1) The original used
% stateSeq = 0 for unlabelled frames. The revised code assigns NOTA = 7
% (one beyond the last labelled state) so the value is a valid positive
% index — this avoids edge cases where 0 could be confused with
% "uninitialised". NOTA frames are never counted in any TPM; they act
% purely as censoring markers that break the transition sequence without
% introducing a spurious 7th state. 
% Partial TPM (labelledTrans) The
% original validTrans mask was fromSeq>0 & toSeq>0, which had the same
% intent but would have included any accidental state-0 frames. The new
% labelledTrans = (fromSeq~=NOTA) & (toSeq~=NOTA) is explicit and robust.
% Row-normalisation in computeTPM now divides by the count of labelled
% departures only, so each row sums to 1 over the labelled states.
% 
% Condition masking interaction The clean way to condition is stateMask =
% stateMask & condition(:) before building stateSeq. Out-of-condition
% frames then become NOTA automatically and are censored by labelledTrans —
% no special interaction with the event-window logic is needed. Event
% windows that straddle a condition boundary simply yield fewer valid
% transitions inside them, which is the correct behaviour. The NOTA rate
% diagnostic printed per event will tell you if a particular event type is
% heavily affected by this. 
% NOTA rate diagnostic Each results(e) struct now
% carries nota_rate_in_window — the fraction of in-window frame-pairs that
% were censored. This is printed in the subplot titles and in a summary
% table. A high NOTA rate for a particular event (say >30%) is a warning
% that your post-event TPM for that event is estimated from sparse labelled
% data and should be interpreted cautiously.

statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
eventnames = {'failed_approach','contact_loss','contact_gain','intercept','cricket_jump','rangemin'};

nStates = numel(statenames);
nEvents = numel(eventnames);
NOTA    = nStates + 1;   % index reserved for none-of-the-above

% ── Build stateMask ─────────────────────────────────────────────────────────
stateMask = logical([hotpursuit(:), chase(:), follow(:), stalk(:), wander(:), pause(:)]);

% Optional condition masking.
% Mask columns of stateMask so out-of-condition frames become NOTA.
% Uncomment one block and set condition / condition_name before running.
%
 c = 3; %set outside in a script, for now
if     c==1; condition = laseron &  dark;  condition_name = 'dark laser on';
elseif c==2; condition = laseron & ~dark;  condition_name = 'light laser on';
elseif c==3; condition = ~laseron &  dark; condition_name = 'dark laser off';
elseif c==4; condition = ~laseron & ~dark; condition_name = 'light laser off';
else
    condition_name = 'all conditions';
end
stateMask = stateMask & condition(:);   % broadcasts condition over columns

% Without a condition, set a neutral label for figure titles:
% condition_name = 'all conditions';

% ── Convert stateMask to categorical sequence ────────────────────────────────
% stateSeq(t) = i    if state i is active at frame t   (1..nStates)
%             = NOTA if no labelled state is active
num_frames = size(stateMask, 1);
stateSeq   = repmat(NOTA, num_frames, 1);
for s = 1:nStates
    stateSeq(stateMask(:,s)) = s;
end
% If stateMask is truly mutually exclusive this loop is fine; the last
% writer wins if there is accidental overlap (add a sanity check below if
% needed).

fprintf('NOTA frames: %d / %d (%.1f%%)\n', ...
    sum(stateSeq==NOTA), num_frames, 100*mean(stateSeq==NOTA));

% ── Event frames ─────────────────────────────────────────────────────────────
fps = 200;
eventFrames = { ...
    failed_approach_event_frames, ...
    contact_loss_event_frames,    ...
    contact_gain_event_frames,    ...
    intercept_event_frames,       ...
    cricket_jump_event_frames,    ...
    rangemin_event_frames         ...
};

% ── Parameters ───────────────────────────────────────────────────────────────
post_event_window_sec = 2;
winFrames = post_event_window_sec*fps;   % post-event window length (e.g. 1 s at 30 fps) -- adjust as needed

% ── Transition pairs ─────────────────────────────────────────────────────────
fromSeq = stateSeq(1:end-1);
toSeq   = stateSeq(2:end);

% Partial TPM: only labelled -> labelled transitions enter counts.
% Transitions touching NOTA are excluded (censored), not treated as a
% transition to/from an extra state.
%
%labelledTrans = (fromSeq ~= NOTA) & (toSeq ~= NOTA); %includes self-transitions (so it shows dwell time)
labelledTrans = (fromSeq ~= NOTA) & (toSeq ~= NOTA) & (fromSeq ~= toSeq); %excludes self-transitions
%labelledTrans =  (fromSeq ~= toSeq); %excludes only self-transitions

% ── Baseline TPM (all valid transitions across all time) ─────────────────────
TPM_baseline = computeTPM(fromSeq(labelledTrans), toSeq(labelledTrans), nStates);

% ── Per-event post-event TPMs ─────────────────────────────────────────────────
results = struct();
for e = 1:nEvents
    % Mark transition-index t as in-window if frame t falls within winFrames
    % of any event of type e.  NOTA frames inside the window are still marked
    % in-window — they just don't contribute counts because labelledTrans
    % will filter them out.
    inWindow = false(num_frames-1, 1);
    for f = eventFrames{e}(:)'
        idx = f : (f + winFrames - 1);
        idx = idx(idx >= 1 & idx <= num_frames-1);
        inWindow(idx) = true;
    end

    postMask = inWindow        & labelledTrans;
    baseMask = (~inWindow)     & labelledTrans;

    TPM_post        = computeTPM(fromSeq(postMask), toSeq(postMask), nStates);
    TPM_base_excl   = computeTPM(fromSeq(baseMask), toSeq(baseMask), nStates);
    counts_post     = computeCounts(fromSeq(postMask), toSeq(postMask), nStates);
    counts_base     = computeCounts(fromSeq(baseMask), toSeq(baseMask), nStates);

    % Per-row G-test and Jensen-Shannon divergence
    pvals = nan(nStates, 1);
    jsdiv = nan(nStates, 1);
    for i = 1:nStates
        obsPost = counts_post(i,:);
        obsBase = counts_base(i,:);
        if sum(obsPost) > 0 && sum(obsBase) > 0
            pvals(i) = gtest_pvalue([obsPost; obsBase]);
            jsdiv(i) = jensenshannon(obsPost/sum(obsPost), obsBase/sum(obsBase));
        end
    end

    % NOTA rate: fraction of in-window frames censored, useful diagnostic
    inWindowFrames  = sum(inWindow);
    notaInWindow    = sum(inWindow & (fromSeq == NOTA | toSeq == NOTA));
    notaRate        = notaInWindow / max(inWindowFrames, 1);

    results(e).event                = eventnames{e};
    results(e).TPM_post             = TPM_post;
    results(e).TPM_baseline_excl    = TPM_base_excl;
    results(e).counts_post          = counts_post;
    results(e).counts_base          = counts_base;
    results(e).pvals                = pvals;
    results(e).jsdiv                = jsdiv;
    results(e).n_post_transitions   = sum(postMask);
    results(e).nota_rate_in_window  = notaRate;
    results(e).condition_name       = condition_name;
end

%% Visualization: difference matrices (TPM_post - TPM_baseline) per event
figure('Position',[100 100 1600 800]);
for e = 1:nEvents
    subplot(2,3,e)
    diffMat = results(e).TPM_post - results(e).TPM_baseline_excl;
    imagesc(diffMat, [-1 1]*max(abs(diffMat(:)),[],'all','omitnan'));
    colormap(gca, redblue());
    colorbar; axis square;

    rownames_star = statenames;
    rownames_star(results(e).pvals < .05) = ...
        append(rownames_star(results(e).pvals < .05), ' *');

    set(gca, 'XTick',1:nStates, 'XTickLabel',statenames, 'XTickLabelRotation',90, ...
             'YTick',1:nStates, 'YTickLabel',rownames_star);
    title(sprintf('%s  (n=%d, NOTA=%.0f%%, win=%d fr)', ...
        eventnames{e}, results(e).n_post_transitions, ...
        100*results(e).nota_rate_in_window, winFrames), 'interp','none');
end
sgtitle(sprintf('TPM_{post} - TPM_{base}  (row=from, col=to)   *p<0.05   [%s]', ...
    condition_name));

%% Summary table
fprintf('\nSignificant row-wise shifts (p<0.05), Jensen-Shannon divergence:\n');
fprintf('%-15s', 'State');
for e = 1:nEvents, fprintf('%-22s', eventnames{e}); end
fprintf('\n');
for i = 1:nStates
    fprintf('%-15s', statenames{i});
    for e = 1:nEvents
        p = results(e).pvals(i);
        j = results(e).jsdiv(i);
        if isnan(p)
            fprintf('%-22s', 'n/a');
        else
            marker = '';
            if p < 0.05, marker = '*'; end
            fprintf('%-22s', sprintf('p=%.3f j=%.4f%s', p, j, marker));
        end
    end
    fprintf('\n');
end

% Print NOTA rates per event as a separate diagnostic
fprintf('\nNOTA censoring rate inside post-event windows:\n');
for e = 1:nEvents
    fprintf('  %-25s %.1f%%\n', eventnames{e}, ...
        100*results(e).nota_rate_in_window);
end

%% Circle diagram
figure('Position',[100 100 1600 1000]);
for e = 1:nEvents
    subplot(2,3,e)
    diffMat = results(e).TPM_post - results(e).TPM_baseline_excl;
    plot_tpm_diff_circle(diffMat, statenames, ...
        'Title', sprintf('%s (n=%d) %s, win=%ds', eventnames{e}, results(e).n_post_transitions, condition_name,post_event_window_sec), ...
        'MinAbsDiff', 0.001);
end
sgtitle(sprintf('TPM_{post} - TPM_{base}  [%s]', ...
    condition_name));

%plot a single figure panel for failed approach
figure
e=1;
diffMat = results(e).TPM_post - results(e).TPM_baseline_excl;
plot_tpm_diff_circle(diffMat, statenames, ...
    'Title', sprintf('TPMpost - TPMbaseline\n%s (n=%d) %s, win=%ds', eventnames{e}, results(e).n_post_transitions, condition_name, post_event_window_sec), ...
    'MinAbsDiff', 0.001);


%% ===== Helper functions =====================================================

function TPM = computeTPM(fromSeq, toSeq, nStates)
% Row-normalise over labelled destinations only (partial TPM).
% Rows with no observed departures are left as zero rather than uniform.
    counts = computeCounts(fromSeq, toSeq, nStates);
    rowsums = sum(counts, 2);
    TPM = counts ./ rowsums;
    TPM(rowsums == 0, :) = 0;
end

function counts = computeCounts(fromSeq, toSeq, nStates)
    counts = zeros(nStates, nStates);
    for k = 1:numel(fromSeq)
        counts(fromSeq(k), toSeq(k)) = counts(fromSeq(k), toSeq(k)) + 1;
    end
end

function p = gtest_pvalue(tbl)
    rowSums  = sum(tbl, 2);
    colSums  = sum(tbl, 1);
    total    = sum(tbl(:));
    expected = rowSums * colSums / total;
    valid    = tbl > 0 & expected > 0;
    G        = 2 * sum(tbl(valid) .* log(tbl(valid) ./ expected(valid)));
    df       = (size(tbl,1)-1) * (size(tbl,2)-1);
    p        = 1 - chi2cdf(G, df);
end

function d = jensenshannon(p, q)
    p = p + eps; q = q + eps;
    p = p/sum(p); q = q/sum(q);
    m = 0.5*(p + q);
    d = sqrt(0.5*sum(p.*log2(p./m)) + 0.5*sum(q.*log2(q./m)));
end

function cmap = redblue()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    cmap = [r' g' b'];
end