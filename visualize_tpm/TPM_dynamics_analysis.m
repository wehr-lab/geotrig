% TPM_dynamics_analysis
%% ============================================================
%  TPM DYNAMICS ANALYSIS
%  Inputs assumed:
%    TPM        : [nStates x nStates] row-stochastic transition matrix
%    statenames : cell array of state name strings
%    fps        : frames per second (for time axis labeling)
% =============================================================
TPM=TPMdark;
% TPM=TPMlight;
fps=200;

nStates = size(TPM, 1);
maxPow  = 200;   % max number of steps to examine

%% ============================================================
%  PART 1A: TPM powers over time
% =============================================================

% Compute TPM^t for t = 1:maxPow
TPMpow = zeros(nStates, nStates, maxPow);
TPMpow(:,:,1) = TPM;
for t = 2:maxPow
    TPMpow(:,:,t) = TPMpow(:,:,t-1) * TPM;
end

% Plot evolution of each entry TPM^t(i,j) vs t
figure('Position',[100 100 1400 900]);
tiledlayout(nStates, nStates, 'TileSpacing','none', 'Padding','compact');
steadyState = squeeze(TPMpow(:,:,end));
for i = 1:nStates
    for j = 1:nStates
        nexttile;
        vals = squeeze(TPMpow(i,j,:));
        plot(1:maxPow, vals, 'k-', 'LineWidth', 0.8); hold on;
        yline(steadyState(i,j), 'r--', 'LineWidth', 0.5);
        ylim([0 1]); xlim([1 maxPow]);
        set(gca,'XTick',[],'YTick',[]);
        if i==1, title(statenames{j},'FontSize',7,'Interpreter','none'); end
        if j==1, ylabel(statenames{i},'FontSize',7,'Interpreter','none','Rotation',45,'HorizontalAlignment','right'); end
    end
end
sgtitle('TPM^t(i,j) vs step t  [row=from, col=to] — red dashed = steady state');

% Also plot as heatmap snapshots at selected time steps
snapshots = round([1, 2, 5, 10, 20, 50, maxPow]);
snapshots = snapshots(snapshots <= maxPow);
figure('Position',[100 100 1600 350]);
tiledlayout(1, numel(snapshots), 'TileSpacing','compact','Padding','compact');
for k = 1:numel(snapshots)
    nexttile;
    imagesc(TPMpow(:,:,snapshots(k)), [0 1]);
    axis square; colorbar;
    set(gca,'XTick',1:nStates,'XTickLabel',statenames,'XTickLabelRotation',90,...
            'YTick',1:nStates,'YTickLabel',statenames,'FontSize',7);
    title(sprintf('t = %d steps (%.1fs)', snapshots(k), snapshots(k)/fps));
end
sgtitle('TPM^t snapshots — rows converging = fast mixing, rows staying different = slow');

%% ============================================================
%  PART 1B: State probability trajectories from each starting state
% =============================================================

figure('Position',[100 100 1400 900]);
cols = lines(nStates);
tiledlayout(ceil(nStates/2), 2, 'TileSpacing','compact','Padding','compact');
for i = 1:nStates
    nexttile;
    hold on;
    for j = 1:nStates
        plot(1:maxPow, squeeze(TPMpow(i,j,:)), '-', 'Color', cols(j,:), 'LineWidth', 1.2);
    end
    title(sprintf('Start: %s', statenames{i}), 'Interpreter','none','FontSize',9);
    xlabel('Steps'); ylabel('P(state at t)');
    ylim([0 1]); xlim([1 maxPow]);
    if i == nStates
        legend(statenames, 'Location','eastoutside','Interpreter','none','FontSize',7);
    end
end
sgtitle('Probability of being in each state at time t, given starting state');

%% ============================================================
%  PART 1C: Relaxation timescales from eigenspectrum
% =============================================================

[V, D] = eig(TPM', 'vector');   % left eigenvectors (columns of V)
[D, sortIdx] = sort(abs(D), 'descend');
V = V(:, sortIdx);

% Relaxation timescales (skip eigenvalue 1 = stationary)
tauRelax = -1 ./ log(D(2:end));
tauRelax(isinf(tauRelax)) = NaN;

figure('Position',[100 100 1000 400]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
stem(1:nStates, D, 'filled', 'LineWidth',1.5);
xlabel('Eigenvalue index'); ylabel('|\lambda|');
xline(1.5,'r--'); yline(1,'k:');
title('Eigenvalue magnitudes'); grid on;
set(gca,'XTick',1:nStates,'XTickLabel',1:nStates);

nexttile;
stem(2:nStates, tauRelax, 'filled', 'LineWidth',1.5);
text((2:nStates)', tauRelax, arrayfun(@(x) sprintf('%.3fs', x/fps), tauRelax, 'UniformOutput', false))
xlabel('Mode index (1 = stationary, skipped)');
ylabel('Relaxation timescale (steps)');
title('Implied timescales  \tau = -1/log|\lambda|'); grid on;
yyaxis right;
ylabel('Time (s)');
ylim([0, max(tauRelax(isfinite(tauRelax)))/fps]);

%% ============================================================
%  PART 1D: Slow eigenvector visualization (state space geometry)
% =============================================================

% Use first 3 non-stationary left eigenvectors (slow modes)
nModes = min(3, nStates-1);
slowVecs = real(V(:, 2:nModes+1));   % skip eigenvec 1 (stationary distribution)

% Stationary distribution (for node sizing)
statDist = abs(V(:,1)); statDist = statDist/sum(statDist);
nodeSize = 60 + 400*statDist;

cols = lines(nStates);
figure('Position',[100 100 1400 500]);
tiledlayout(1, nModes-1+1, 'TileSpacing','compact','Padding','compact');

% 2D projection: mode 2 vs mode 3
if nModes >= 2
    nexttile;
    hold on;
    for i = 1:nStates
        scatter(slowVecs(i,1), slowVecs(i,2), nodeSize(i), cols(i,:), 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8);
        text(slowVecs(i,1), slowVecs(i,2), sprintf('  %s',statenames{i}), ...
            'FontSize',9,'Interpreter','none');
    end
    xlabel('Slow mode 1 (\lambda_2)'); ylabel('Slow mode 2 (\lambda_3)');
    title('State geometry: slow eigenvectors'); grid on; axis equal;
end

% 3D projection: modes 2,3,4
if nModes >= 3
    nexttile;
    hold on;
    for i = 1:nStates
        scatter3(slowVecs(i,1), slowVecs(i,2), slowVecs(i,3), nodeSize(i), ...
            cols(i,:), 'filled','MarkerEdgeColor','k','LineWidth',0.8);
        text(slowVecs(i,1), slowVecs(i,2), slowVecs(i,3), ...
            sprintf('  %s',statenames{i}), 'FontSize',9,'Interpreter','none');
    end
    xlabel('Mode 1'); ylabel('Mode 2'); zlabel('Mode 3');
    title('State geometry: 3 slow modes'); grid on; view(45,30);
end

%% ============================================================
%  PART 2: METASTABILITY — Almost-invariant sets (PCCA-inspired)
% =============================================================
% Uses k-means on slow eigenvectors to find metastable clusters
% (full PCCA+ requires specialized toolbox; this is a practical approximation)

nClusters = 2;   % adjust based on spectral gap — look for gap in eigenvalue plot

% Normalize rows of slow eigenvector matrix before clustering (Shi-Malik style)
X = slowVecs(:, 1:min(nModes, nStates-1));
X = X ./ (sqrt(sum(X.^2, 2)) + eps);  % row-normalize

% k-means clustering in slow eigenvector space
rng(42);
[clusterIdx, clusterCenters] = kmeans(X, nClusters, 'Replicates', 20, 'Distance','cosine');

% Compute coarse-grained TPM over metastable sets
TPM_coarse = zeros(nClusters, nClusters);
for a = 1:nClusters
    for b = 1:nClusters
        statesA = find(clusterIdx == a);
        statesB = find(clusterIdx == b);
        TPM_coarse(a,b) = sum(sum(statDist(statesA) .* TPM(statesA, statesB)));
    end
    if sum(TPM_coarse(a,:)) > 0
        TPM_coarse(a,:) = TPM_coarse(a,:) / sum(TPM_coarse(a,:));
    end
end

% Metastability index: mean of diagonal (probability of staying in same macrostate)
metastability = trace(TPM_coarse) / nClusters;
fprintf('\nMetastability index (%.0f clusters): %.3f  (1=perfectly metastable, 1/n=random)\n', ...
    nClusters, metastability);

% Print cluster membership
fprintf('\nMetastable cluster assignments:\n');
for c = 1:nClusters
    members = statenames(clusterIdx == c);
    fprintf('  Cluster %d: %s\n', c, strjoin(members, ', '));
end

% Visualize clusters on slow eigenvector plot
clusterCols = [0.85 0.2 0.2; 0.2 0.6 0.85; 0.3 0.75 0.4; 0.9 0.6 0.1; 0.6 0.3 0.8];
figure('Position',[100 100 1400 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
for i = 1:nStates
    c = clusterIdx(i);
    scatter(slowVecs(i,1), slowVecs(i,2), nodeSize(i), clusterCols(c,:), ...
        'filled','MarkerEdgeColor','k','LineWidth',1.2);
    text(slowVecs(i,1), slowVecs(i,2), sprintf('  %s',statenames{i}), ...
        'FontSize',9,'Interpreter','none');
end
xlabel('Slow mode 1'); ylabel('Slow mode 2');
title('Metastable clusters in eigenvector space'); grid on; axis equal;

nexttile;
clusterLabels = arrayfun(@(c) sprintf('Cluster %d',c), 1:nClusters, 'UniformOutput',false);
imagesc(TPM_coarse, [0 1]); colorbar; axis square;
set(gca,'XTick',1:nClusters,'XTickLabel',clusterLabels,'XTickLabelRotation',45,...
        'YTick',1:nClusters,'YTickLabel',clusterLabels);
title(sprintf('Coarse-grained TPM over metastable sets\nMetastability index = %.3f', metastability));
colormap(gca, hot);
for a = 1:nClusters
    for b = 1:nClusters
        text(b, a, sprintf('%.2f',TPM_coarse(a,b)), 'HorizontalAlignment','center',...
            'FontSize',10,'Color', TPM_coarse(a,b)>0.5*[1 1 1]+(TPM_coarse(a,b)<=0.5)*[0 0 0]);
    end
end

str={...
'**Key things to look at in the output:**',...
'',...
'- **Eigenvalue plot (1C):** the spectral gap (big drop between λ₂ and λ₃, or λ₃ and λ₄) tells you how many metastable clusters are meaningful — set `nClusters` to match.',...
'- **Implied timescales (1C right):** if τ₂ >> τ₃, there is one dominant slow process; if τ₂ ≈ τ₃ >> τ₄, there are two comparable slow modes.',...
'- **Eigenvector scatter (1D):** states that cluster tightly transition rapidly among themselves; states far apart mix slowly — this is the geometric picture of your dynamics.',...
'- **Metastability index (2):** ranges from 1/nClusters (random, no metastability) to 1 (perfectly trapped). Values above ~0.7 indicate genuine metastable structure worth interpreting.'};

% conclusions:
% there's one dominant slow process
% there are 2 clusters, so I set nClusters accordingly

% chase1 and chase2 cluster tightly, thus they transition rapidly among themselves;
% pause and stalk cluster tightly, thus they transition rapidly among themselves;
% states far apart mix slowly —
% wander, and chase3 transition slowly with each other and the rest
% chase3 slowly dominates
% this is the geometric picture of your dynamics.

% Metastability index (2 clusters): 0.500
% this is lower than .7 
% Values above ~0.7 indicate genuine metastable structure worth interpreting.





