function plot_tpm_sankey(P, nSteps, varargin)
% PLOT_TPM_SANKEY  Visualize multi-step probability flow under a TPM as a
% Sankey-style diagram: one column of nodes per time step, with flows
% between columns proportional to the probability mass moving from state
% i (at time t) to state j (at time t+1).
%
% USAGE:
%   plot_tpm_sankey(P, nSteps)
%   plot_tpm_sankey(P, nSteps, 'pi0', pi0)
%   plot_tpm_sankey(P, nSteps, 'partition', partition)
%
% INPUTS:
%   P      - N x N transition probability matrix (rows sum to 1)
%   nSteps - number of time steps to show (number of columns = nSteps+1)
%
% OPTIONAL NAME-VALUE ARGS:
%   'pi0'         - 1 x N initial distribution (default: point mass on
%                   state 1)
%   'partition'   - N x 1 vector of cluster labels (1..K). If provided,
%                   nodes are grouped vertically by cluster (with small
%                   gaps between groups) and colored by cluster, using
%                   the SAME [~,~,partIdx]=unique(partition) convention
%                   as the other functions. Flows are still drawn at the
%                   fine-state level, but coloring makes within-cluster
%                   vs. cross-cluster flows easy to distinguish.
%   'stateNames'  - cell array of N strings for fine-state labels
%   'minFlow'     - flows (in probability mass) below this are not drawn
%                   (default 0.01)
%   'colors'      - N x 3 matrix of RGB colors for fine states (default:
%                   if partition given, derived from cluster colors via
%                   lines(K); otherwise lines(N))
%
% NOTES:
%   - Node height at column t represents pi_t(state) = (pi0 * P^t)(state).
%   - Flow (i,t) -> (j,t+1) has "thickness" pi_t(i) * P(i,j), i.e. the
%     probability mass transitioning from i to j between steps t and t+1.
%   - This complements plot_occupancy_over_time: occupancy shows total
%     mass per state over time, while this shows HOW that mass
%     redistributes step by step.
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];
%   plot_tpm_sankey(P, 5, 'partition', partition, ...
%       'stateNames', {'A1','A2','B1','B2'}, 'pi0', [1 0 0 0]);

    p = inputParser;
    addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addRequired(p, 'nSteps', @(x) isscalar(x) && x>=1);
    addParameter(p, 'pi0', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'partition', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'stateNames', {}, @iscell);
    addParameter(p, 'minFlow', 0.01, @(x) isscalar(x) && x>=0);
    addParameter(p, 'colors', [], @(x) isempty(x) || (ismatrix(x) && size(x,2)==3));
    parse(p, P, nSteps, varargin{:});

    P = p.Results.P;
    nSteps = p.Results.nSteps;
    N = size(P,1);

    if isempty(p.Results.stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    else
        stateNames = p.Results.stateNames;
    end

    if isempty(p.Results.pi0)
        pi0 = zeros(1,N);
        pi0(1) = 1;
    else
        pi0 = p.Results.pi0(:)';
        if abs(sum(pi0)-1) > 1e-6
            warning('plot_tpm_sankey:pi0NotNormalized', ...
                'pi0 does not sum to 1 (sum = %.6f); normalizing.', sum(pi0));
            pi0 = pi0/sum(pi0);
        end
    end

    minFlow = p.Results.minFlow;
    partition = p.Results.partition;

    % --- Determine vertical ordering & colors ---
    if isempty(partition)
        order = 1:N;
        partIdx = ones(N,1);
        K = 1;
        gapEach = 0;
    else
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_tpm_sankey:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);
        K = max(partIdx);
        [~, order] = sort(partIdx); % cluster-major ordering
        gapEach = 0.5; % vertical gap (in "probability units") between clusters
    end

    if isempty(p.Results.colors)
        if K > 1
            clusterColors = lines(K);
            colors = clusterColors(partIdx, :); % N x 3, one row per fine state
        else
            colors = lines(N);
        end
    else
        colors = p.Results.colors;
    end

    % --- Propagate distribution forward ---
    pi_history = zeros(nSteps+1, N);
    pi_history(1,:) = pi0;
    for t = 1:nSteps
        pi_history(t+1,:) = pi_history(t,:) * P;
    end

    % --- Compute vertical layout (y-positions for each state at each column) ---
    % Each column's total height = 1 (total probability) + gaps between clusters.
    % We compute, for each state in `order`, its [yTop, yBottom] band at each
    % time step, stacked top-to-bottom according to `order`, with small gaps
    % inserted between clusters.

    % Precompute cumulative gap offsets for cluster boundaries
    clusterOfOrdered = partIdx(order); % K-grouped order
    gapBefore = zeros(N,1); % extra vertical offset before this state (cumulative)
    cumGap = 0;
    for ii = 1:N
        if ii > 1 && clusterOfOrdered(ii) ~= clusterOfOrdered(ii-1)
            cumGap = cumGap + gapEach;
        end
        gapBefore(ii) = cumGap;
    end
    totalHeight = 1 + cumGap;

    % yTopOrdered(t, ii) = top y-coordinate of the ii-th state (in `order`)
    % at time step t (columns indexed 1..nSteps+1)
    yTop = zeros(nSteps+1, N);
    yBot = zeros(nSteps+1, N);
    for t = 1:nSteps+1
        cursor = totalHeight; % start from top, go downward
        for ii = 1:N
            i = order(ii);
            h = pi_history(t,i);
            yTop(t,ii) = cursor - gapBefore(ii)*0 ; % placeholder, adjusted below
        end
    end

    % Recompute properly: cursor decreases by node height + gap-if-new-cluster
    for t = 1:nSteps+1
        cursor = totalHeight;
        prevCluster = -1;
        for ii = 1:N
            i = order(ii);
            cl = clusterOfOrdered(ii);
            if prevCluster ~= -1 && cl ~= prevCluster
                cursor = cursor - gapEach;
            end
            h = pi_history(t,i);
            yTop(t,ii) = cursor;
            yBot(t,ii) = cursor - h;
            cursor = cursor - h;
            prevCluster = cl;
        end
    end

    % --- Plot ---
    figure('Color','w');
    hold on; axis off;

    xCols = 1:(nSteps+1);
    nodeWidth = 0.15;

    % --- Draw flows first (so nodes sit on top) ---
    for t = 1:nSteps
        % For each pair (i ordered as ii, j ordered as jj), flow mass =
        % pi_history(t,i) * P(i,j)
        % We need to track, within each source node's band, how much of
        % its height has already been allocated to outgoing flows (so
        % multiple flows from the same node stack without overlapping),
        % and similarly for incoming flows into each destination node.

        outCursor = yTop(t,:);   % current top of remaining unallocated band, per source (ordered index)
        inCursor  = yTop(t+1,:); % current top of remaining unallocated band, per dest (ordered index)

        for ii = 1:N
            i = order(ii);
            if pi_history(t,i) < eps
                continue;
            end
            for jj = 1:N
                j = order(jj);
                flow = pi_history(t,i) * P(i,j);
                if flow < minFlow
                    continue;
                end

                ySrcTop = outCursor(ii);
                ySrcBot = ySrcTop - flow;
                outCursor(ii) = ySrcBot;

                yDstTop = inCursor(jj);
                yDstBot = yDstTop - flow;
                inCursor(jj) = yDstBot;

                draw_flow_ribbon(xCols(t)+nodeWidth/2, xCols(t+1)-nodeWidth/2, ...
                    ySrcTop, ySrcBot, yDstTop, yDstBot, colors(i,:));
            end
        end
    end

    % --- Draw nodes ---
    for t = 1:nSteps+1
        for ii = 1:N
            i = order(ii);
            h = pi_history(t,i);
            if h < eps
                continue;
            end
            x0 = xCols(t) - nodeWidth/2;
            rectangle('Position', [x0, yBot(t,ii), nodeWidth, h], ...
                'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 0.5);

            % Label state name next to node (only at first and last columns
            % to avoid clutter)
            if t == 1
                text(x0 - 0.03, (yTop(t,ii)+yBot(t,ii))/2, stateNames{i}, ...
                    'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
                    'FontSize', 9);
            end
            if t == nSteps+1
                text(x0 + nodeWidth + 0.03, (yTop(t,ii)+yBot(t,ii))/2, stateNames{i}, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                    'FontSize', 9);
            end
        end
    end

    % --- Time step labels along the bottom ---
    for t = 1:nSteps+1
        text(xCols(t), -0.05, sprintf('t=%d', t-1), ...
            'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 10);
    end

    xlim([0.5, nSteps+1.5]);
    ylim([-0.15, totalHeight+0.05]);
    title('Multi-Step Probability Flow (Sankey Diagram)');
    hold off;
end

% ===================== HELPER FUNCTIONS =====================

function draw_flow_ribbon(x1, x2, y1top, y1bot, y2top, y2bot, color)
% Draws a smooth (sigmoid-interpolated) ribbon from a vertical band
% [y1bot,y1top] at x=x1 to a vertical band [y2bot,y2top] at x=x2.

    nPts = 30;
    xs = linspace(x1, x2, nPts);

    % Smooth interpolation (cosine/sigmoid-like) for top and bottom edges
    s = (1 - cos(linspace(0,pi,nPts))) / 2; % 0 -> 1 smoothly
    topEdge = y1top + (y2top - y1top) * s;
    botEdge = y1bot + (y2bot - y1bot) * s;

    xPoly = [xs, fliplr(xs)];
    yPoly = [topEdge, fliplr(botEdge)];

    fill(xPoly, yPoly, color, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
end
