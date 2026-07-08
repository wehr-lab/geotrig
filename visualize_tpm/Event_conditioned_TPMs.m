% Event_conditioned_TPMs

% maybe replace rangemin with intercept, or at least rename it
% add escape?

% 1. Event-conditioned TPMs (stratification)
% 
% Split transitions into two pools — those occurring within some window
% after an event, vs. baseline — and build separate TPMs for each. Compare
% with a chi-square or G-test on the transition count tables, or row-by-row
% comparison of post-event vs. baseline distributions (e.g. KL divergence
% or Jensen-Shannon distance per row). Simple, interpretable, and a good
% first pass to establish whether events matter before building anything
% fancier.

statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
eventnames = {'failed_approach','contact_loss','contact_gain','intercept','cricket_jump', 'rangemin'};

nStates = numel(statenames);
nEvents = numel(eventnames);

% --- Inputs expected ---
% stateMask : [num_frames x nStates] logical, mutually exclusive (post-refinement)
% stateMask(t,i) = true if state i is active at time t
% eventFrames : cell array, eventFrames{e} = list of frame indices for event e
% fps : frame rate (for converting window to frames), or just set winFrames directly

stateMask = [    hotpursuit(:),chase(:),follow(:), stalk(:), wander(:), pause(:)];
%stateMask = [  chase(:), stalk(:), wander(:), pause(:)];

c=3;
    if c==1
        condition=laseron & dark; condition_name{c}='dark laser on';
    elseif c==2
        condition=laseron & light; condition_name{c}='light laser on';
    elseif c==3
        condition=~laseron & dark; condition_name{c}='dark laser off';
    elseif c==4
        condition=~laseron & light; condition_name{c}='light laser off';
    end


stateMask = [    condition & hotpursuit(:), condition & chase(:), condition & follow(:),...
     condition & stalk(:),  condition & wander(:),  condition & pause(:)];


size(stateMask)
fps=200;

eventFrames={...
    failed_approach_event_frames,...
    contact_loss_event_frames,...
    contact_gain_event_frames,...
    intercept_event_frames, ...
    cricket_jump_event_frames,...
    rangemin_event_frames,...
    };

% --- Parameters ---
post_event_window_sec = 2;
winFrames = post_event_window_sec*fps;   % post-event window length (e.g. 1 s at 30 fps) -- adjust as needed

num_frames = size(stateMask,1);

% Convert state mask to a single categorical sequence (1..nStates), 0 = none active
stateSeq = zeros(num_frames,1);
[hasState, stateSeq(stateSeq==0)] = max(stateMask,[],2);
stateSeq(~any(stateMask,2)) = 0;

%% Build baseline TPM (all transitions) and per-event post-event TPMs
% Transitions occur between frame t and t+1
fromSeq = stateSeq(1:end-1);
toSeq   = stateSeq(2:end);
% validTrans = fromSeq>0 & toSeq>0;  % exclude frames where no state active
%note that this excludes transitions to/from "none of the above"
validTrans = fromSeq>0 & toSeq>0 & (fromSeq ~= toSeq);  % exclude frames where no state active AND exclude self-transitions

% Baseline TPM (all valid transitions)
TPM_baseline = computeTPM(fromSeq(validTrans), toSeq(validTrans), nStates);

% Mark which transition indices (t -> t+1, i.e., index t) fall in a post-event window
results = struct();
for e = 1:nEvents
    inWindow = false(num_frames-1,1);
    for f = eventFrames{e}(:)'
        idx = f:(f+winFrames-1);
        idx = idx(idx>=1 & idx<=num_frames-1);
        inWindow(idx) = true;
    end

    postMask = inWindow & validTrans;
    baseMask = ~inWindow & validTrans;

    TPM_post = computeTPM(fromSeq(postMask), toSeq(postMask), nStates);
    TPM_base_excl = computeTPM(fromSeq(baseMask), toSeq(baseMask), nStates);

    % Counts for stats
    counts_post = computeCounts(fromSeq(postMask), toSeq(postMask), nStates);
    counts_base = computeCounts(fromSeq(baseMask), toSeq(baseMask), nStates);

    % Per-row chi-square / G-test comparing post vs baseline distribution
    pvals = nan(nStates,1);
    jsdiv = nan(nStates,1);
    for i = 1:nStates
        obsPost = counts_post(i,:);
        obsBase = counts_base(i,:);
        if sum(obsPost)>0 && sum(obsBase)>0
            % G-test (log-likelihood ratio) via contingency table
            tbl = [obsPost; obsBase];
            pvals(i) = gtest_pvalue(tbl);
            jsdiv(i) = jensenshannon(obsPost/sum(obsPost), obsBase/sum(obsBase));
        end
    end

    results(e).event = eventnames{e};
    results(e).TPM_post = TPM_post;
    results(e).TPM_baseline_excl = TPM_base_excl;
    results(e).counts_post = counts_post;
    results(e).counts_base = counts_base;
    results(e).pvals = pvals;
    results(e).jsdiv = jsdiv;
    results(e).n_post_transitions = sum(postMask);
end

%% Visualization: difference matrices (TPM_post - TPM_baseline) per event
figure('Position',[100 100 1600 800]);
for e = 1:nEvents
    subplot(2,3,e)
    diffMat = results(e).TPM_post - results(e).TPM_baseline_excl;
    imagesc(diffMat, [-1 1]*max(abs(diffMat(:)),[],'all','omitnan'));
    colormap(gca, redblue());
    colorbar; axis square;

    % add stars to diff matrices
    rownames_star=statenames;
    rownames_star(results(e).pvals<.05)=append(rownames_star(results(e).pvals<.05), ' *');

    set(gca,'XTick',1:nStates,'XTickLabel',statenames,'XTickLabelRotation',90, ...
        'YTick',1:nStates,'YTickLabel',rownames_star);
    title(sprintf('%s (n=%d, win=%d fr)', eventnames{e}, results(e).n_post_transitions, winFrames), 'interp', 'none');
end
sgtitle(['TPM_{post-event} - TPM_{baseline}  (row = from-state, col = to-state), *=significant row-wise shifts (p<0.05)', condition_name{c}]);



%% Summary table: which rows are significantly affected by which events
fprintf('\nSignificant row-wise shifts (p<0.05), with Jensen-Shannon divergence:\n');
fprintf('%-15s', 'State');
for e = 1:nEvents, fprintf('%-18s', eventnames{e}); end
fprintf('\n');
for i = 1:nStates
    fprintf('%-15s', statenames{i});
    for e = 1:nEvents
        p = results(e).pvals(i);
        j = results(e).jsdiv(i);
        if isnan(p)
            fprintf('%-18s','n/a');
        else
            marker = '';
            if p < 0.05, marker = '*'; end
            fprintf('%-18s', sprintf('p=%.3f j=%.4f%s', p, j, marker));
        end
    end
    fprintf('\n');
end



%% circle diagram of difference TPMs
figure('Position',[100 100 1600 1000]);
for e = 1:nEvents
    subplot(2,3,e)
    diffMat = results(e).TPM_post - results(e).TPM_baseline_excl;
    plot_tpm_diff_circle(diffMat, statenames, ...
        'Title', sprintf('%s (n=%d) %s, win=%ds', eventnames{e}, results(e).n_post_transitions, condition_name{c},post_event_window_sec), ...
        'MinAbsDiff', 0.001);
end


%% ===== Helper functions =====
function TPM = computeTPM(fromSeq, toSeq, nStates)
    counts = computeCounts(fromSeq, toSeq, nStates);
    TPM = counts ./ sum(counts,2);
    TPM(isnan(TPM)) = 0;
end

function counts = computeCounts(fromSeq, toSeq, nStates)
    counts = zeros(nStates, nStates);
    for k = 1:numel(fromSeq)
        counts(fromSeq(k), toSeq(k)) = counts(fromSeq(k), toSeq(k)) + 1;
    end
end

function p = gtest_pvalue(tbl)
    % tbl: 2 x nCategories contingency table
    rowSums = sum(tbl,2);
    colSums = sum(tbl,1);
    total = sum(tbl(:));
    expected = rowSums * colSums / total;
    valid = tbl>0 & expected>0;
    G = 2 * sum(tbl(valid) .* log(tbl(valid)./expected(valid)));
    df = (size(tbl,1)-1)*(size(tbl,2)-1);
    p = 1 - chi2cdf(G, df);
end

function d = jensenshannon(p, q)
    p = p + eps; q = q + eps;
    p = p/sum(p); q = q/sum(q);
    m = 0.5*(p+q);
    d = sqrt(0.5*sum(p.*log2(p./m)) + 0.5*sum(q.*log2(q./m)));
end

function cmap = redblue()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    cmap = [r' g' b'];
end



% Notes / things to adjust:
% 
% winFrames: set based on your frame rate and expected effect timescale
% (this becomes the basis for step 2's time-since-event analysis).
% stateSeq==0 frames (no state active, if your refined states still have
% rare gaps) are excluded from transition counts via validTrans. If you've
% truly eliminated overlap/gaps, this should rarely trigger. The G-test
% (gtest_pvalue) compares each row's post-event transition distribution
% against its baseline distribution — significant p-values flag states
% whose outgoing transition probabilities shift after a given event.
% Jensen-Shannon divergence gives an effect-size companion to the p-value
% (bounded 0-1, robust to zero counts). redblue() is a minimal diverging
% colormap; swap for colormap(redbluecmap) if you have the Bioinformatics
% Toolbox, or colormap(crameri('vik')) etc. if available. Overlapping event
% windows (two events close together) will cause some post-event frames to
% be counted for both events — that's fine for per-event analysis but worth
% checking n_post_transitions isn't dominated by overlap if events cluster.
% 
% This gives you, per event: a difference matrix showing which transitions
% change, plus a significance/effect-size summary identifying which
% from-states are most perturbed by which events — directly motivating
% which states/events to focus on in step (2).
% 
% 
