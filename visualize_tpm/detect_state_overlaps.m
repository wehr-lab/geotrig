%detect state overlaps
%
% Our states are not guaranteed to be exclusive or exhaustive, which
% fundamentally invalidates the TPM. If overlaps are rare, it might not be
% that much of a problem. But approach/chase overlap 100%.
%
% how much overlap is there among all the states?



%% Here's code to quantify and visualize overlap, assuming each state is
% represented as a logical/binary vector over time (1 = state active at
% that sample).

statenames = { ...
     'chase1', 'chase2', 'chase3', 'chase4',  ...
    'stalk', 'pause', 'wander'};

nStates = numel(statenames);

% stateMask: [nSamples x nStates] logical matrix, one column per state
% stateMask(t,i) = true if state i is active at time t

stateMask = [ ...
    chase1(:),chase2(:),chase3(:),chase4(:),...
    stalk(:), pause(:), wander(:)];
size(stateMask)


%% 1. Pairwise overlap quantification
% Jaccard index: |A∩B| / |A∪B|
jaccard = zeros(nStates);
% Overlap coefficient: |A∩B| / min(|A|,|B|)  (fraction of smaller state covered)
overlapCoef = zeros(nStates);
% Conditional: P(B active | A active)
condProb = zeros(nStates);

for i = 1:nStates
    Ai = stateMask(:,i);
    for j = 1:nStates
        Aj = stateMask(:,j);
        inter = sum(Ai & Aj);
        uni   = sum(Ai | Aj);
        jaccard(i,j) = inter / max(uni, eps);
        overlapCoef(i,j) = inter / max(min(sum(Ai), sum(Aj)), eps);
        condProb(i,j) = inter / max(sum(Ai), eps);  % P(j | i)
    end
end

%% 2. Summary stats
nActivePerSample = sum(stateMask, 2);
fprintf('frames with 0 states active: %.2f%%\n', 100*mean(nActivePerSample==0));
fprintf('frames with 1 state active:  %.2f%%\n', 100*mean(nActivePerSample==1));
fprintf('frames with >1 state active: %.2f%%\n', 100*mean(nActivePerSample>1));
fprintf('Max simultaneous states: %d\n', max(nActivePerSample));

%% 3. Heatmaps
figure('Position',[100 100 1400 450]);

subplot(1,3,1)
imagesc(jaccard); axis square; colorbar; clim([0 1]);
set(gca,'XTick',1:nStates,'XTickLabel',statenames,'XTickLabelRotation',90, ...
    'YTick',1:nStates,'YTickLabel',statenames);
title('Jaccard Index (|A∩B|/|A∪B|)');

subplot(1,3,2)
imagesc(overlapCoef); axis square; colorbar; clim([0 1]);
set(gca,'XTick',1:nStates,'XTickLabel',statenames,'XTickLabelRotation',90, ...
    'YTick',1:nStates,'YTickLabel',statenames);
title('Overlap Coefficient (|A∩B|/min(|A|,|B|))');

subplot(1,3,3)
imagesc(condProb); axis square; colorbar; clim([0 1]);
set(gca,'XTick',1:nStates,'XTickLabel',statenames,'XTickLabelRotation',90, ...
    'YTick',1:nStates,'YTickLabel',statenames);
title('P(column state | row state)');

%% 4. Histogram of co-active state counts
figure;
histogram(nActivePerSample, 'BinMethod','integers');
xlabel('Number of simultaneously active states');
ylabel('frame count');
title('Distribution of state co-activation');
xticks(0: max(nActivePerSample));

%% 5. (Optional) UpSet-style table of common combinations
combos = unique(stateMask, 'rows');
comboCounts = zeros(size(combos,1),1);
for k = 1:size(combos,1)
    comboCounts(k) = sum(all(stateMask == combos(k,:), 2));
end
[sortedCounts, idx] = sort(comboCounts, 'descend');
% topN = min(10, numel(idx));
fprintf('\nlist of ranked state combinations:\n');
for k = 1:numel(idx)
    activeIdx = find(combos(idx(k),:));
    if isempty(activeIdx)
        label = '(none active)';
        %skip
    elseif length(activeIdx)==1
        %skip
                label = strjoin(statenames(activeIdx), ' + ');
    elseif length(activeIdx)>1
        label = strjoin(statenames(activeIdx), ' + ');
        fprintf('%6d frames (%.2f%%): %s\n', sortedCounts(k), 100*sortedCounts(k)/numel(nActivePerSample), label);
    else
        fprintf('wtf')
    end
end


















