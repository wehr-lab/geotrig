function [MFPT, fig] = plot_mfpt_heatmap(P, stateNames, varargin)
% PLOT_MFPT_HEATMAP  Compute the mean first-passage time (MFPT) matrix for
% a TPM and display it as a heatmap, optionally reordered/grouped by a
% candidate partition. Block-diagonal structure (short MFPT within
% blocks, long MFPT between blocks) after reordering is a visual
% indicator that the partition aligns with the chain's dynamic structure.
%
% USAGE:
%   MFPT = plot_mfpt_heatmap(P)
%   MFPT = plot_mfpt_heatmap(P, stateNames)
%   MFPT = plot_mfpt_heatmap(P, stateNames, 'partition', partition)
%
% INPUTS:
%   P          - N x N transition probability matrix (rows sum to 1)
%   stateNames - cell array of N strings (use {} for default S1, S2, ...)
%
% OPTIONAL NAME-VALUE ARGS:
%   'partition' - N x 1 vector of cluster labels (1..K). If provided,
%                 states are REORDERED (grouped by cluster) before
%                 plotting, with white grid lines marking cluster
%                 boundaries, so block structure is easy to see.
%
% OUTPUT:
%   MFPT - N x N matrix, MFPT(i,j) = expected number of steps to first
%          reach state j starting from state i (MFPT(i,i) = 0 by
%          convention here)
%   fig  - handle to the created figure
%
% METHOD:
%   Computed via the standard linear system: for each target j, solve
%       m_i = 1 + sum_{k != j} P(i,k) * m_k    for all i != j
%       m_j = 0
%   This is solved column-by-column (one linear solve per target state).
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];
%   MFPT = plot_mfpt_heatmap(P, {'A1','A2','B1','B2'}, 'partition', partition);

    p = inputParser;
    addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addOptional(p, 'stateNames', {}, @iscell);
    addParameter(p, 'partition', [], @(x) isempty(x) || isvector(x));
    parse(p, P, stateNames, varargin{:});

    P = p.Results.P;
    N = size(P,1);

    if isempty(p.Results.stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    else
        stateNames = p.Results.stateNames;
    end

    % --- Compute MFPT matrix ---
    MFPT = zeros(N,N);
    I = eye(N);

    for j = 1:N
        keep = setdiff(1:N, j);
        Psub = P(keep, keep);
        A = I(keep,keep) - Psub;
        b = ones(numel(keep),1);
        m = A \ b; % solve (I - P_sub) m = 1
        MFPT(keep, j) = m;
        MFPT(j,j) = 0;
    end

    % --- Determine display order ---
    partition = p.Results.partition;
    if isempty(partition)
        order = 1:N;
        boundaries = [];
    else
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_mfpt_heatmap:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);
        [sortedIdx, order] = sort(partIdx);

        % Find boundaries between clusters in the reordered matrix
        boundaries = find(diff(sortedIdx) ~= 0);
    end

    MFPT_display = MFPT(order, order);
    names_display = stateNames(order);

    % --- Plot heatmap ---
    fig = figure('Color','w');
    imagesc(MFPT_display);
    colormap(flipud(parula)); % darker = shorter MFPT (more "connected")
    cb = colorbar;
    ylabel(cb, 'Mean first-passage time (steps)');

    set(gca, 'XTick', 1:N, 'XTickLabel', names_display, 'XTickLabelRotation', 45, ...
             'YTick', 1:N, 'YTickLabel', names_display);
    xlabel('To state j');
    ylabel('From state i');
    title('Mean First-Passage Time Matrix');
    axis square;

    % --- Overlay text values ---
    hold on;
    for i = 1:N
        for j = 1:N
            if i == j
                continue;
            end
            text(j, i, sprintf('%.1f', MFPT_display(i,j)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize', 8, 'Color', [0.9 0.9 0.9]);
        end
    end

    % --- Draw cluster boundary lines ---
    if ~isempty(boundaries)
        for b = boundaries(:)'
            xline_pos = b + 0.5;
            plot([xline_pos xline_pos], [0.5 N+0.5], 'w-', 'LineWidth', 2.5);
            plot([0.5 N+0.5], [xline_pos xline_pos], 'w-', 'LineWidth', 2.5);
        end
    end
    hold off;
end
