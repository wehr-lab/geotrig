function plot_tpm_circle_grouped(P, stateNames, partition, opts)
% PLOT_TPM_CIRCLE_GROUPED  Visualize a TPM as a circular state diagram,
% with arrow thickness/darkness encoding transition probability (as in
% plot_tpm_circle), but with states arranged in angular clusters
% according to a partition. This makes it easy to visually assess whether
% a candidate lumping "carves nature at its joints": within-cluster
% arrows should be dark/thick, between-cluster arrows light/thin, with
% clear angular gaps separating the clusters.
%
% USAGE:
%   plot_tpm_circle_grouped(P, stateNames, partition)
%   plot_tpm_circle_grouped(P, stateNames, partition, opts)
%   plot_tpm_circle_grouped(P, stateNames, [], opts)   % no grouping (falls
%                                                       % back to evenly
%                                                       % spaced layout)
%
% INPUTS:
%   P          - N x N transition probability matrix (rows sum to 1)
%   stateNames - cell array of N strings (use {} for default S1, S2, ...)
%   partition  - N x 1 vector of cluster labels (1..K). Pass [] for no
%                grouping (equivalent to plot_tpm_circle's layout).
%   opts       - (optional) struct, same fields as plot_tpm_circle, plus:
%       .clusterGapFrac - fraction of the circle (0..1) devoted to gaps
%                         between clusters, total across all gaps
%                         (default 0.15). Larger = more visual separation.
%       .clusterNames   - cell array of K strings, labeled at the outer
%                         arc of each cluster group (default none)
%       .clusterColors  - K x 3 matrix of RGB colors for cluster arc
%                         labels/backgrounds (default: lines(K))
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];
%   plot_tpm_circle_grouped(P, {'A1','A2','B1','B2'}, partition, ...
%       struct('clusterNames', {{'Cluster A','Cluster B'}}));

    N = size(P,1);

    if nargin < 2 || isempty(stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    end
    if nargin < 3
        partition = [];
    end
    if nargin < 4
        opts = struct();
    end
    if ~isfield(opts,'minProb'),     opts.minProb     = 0.01; end
    if ~isfield(opts,'radius'),      opts.radius      = 1;    end
    if ~isfield(opts,'nodeSize'),    opts.nodeSize    = 60;   end
    if ~isfield(opts,'colormapFn'),  opts.colormapFn  = @(x) (1-x)*[1 1 1]; end
    if ~isfield(opts,'selfLoop'),    opts.selfLoop    = true; end
    if ~isfield(opts,'clusterGapFrac'), opts.clusterGapFrac = 0.15; end
    if ~isfield(opts,'clusterNames'), opts.clusterNames = {}; end
    if ~isfield(opts,'clusterColors'), opts.clusterColors = []; end

    % --- Determine angular position of each state ---
    if isempty(partition)
        % Fall back to evenly spaced layout (no grouping)
        theta_states = compute_even_angles(N);
        partIdx = ones(N,1);
        K = 1;
        groupOrder = 1:N; % identity ordering
    else
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_tpm_circle_grouped:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);
        K = max(partIdx);

        [theta_states, groupOrder, groupAngleRanges] = compute_grouped_angles(partIdx, opts.clusterGapFrac);
    end

    x = opts.radius * cos(theta_states);
    y = opts.radius * sin(theta_states);

    figure('Color','w');
    hold on; axis equal off;

    % --- Optional: draw cluster background wedges/arcs behind everything ---
    if ~isempty(partition) && K > 1
        if isempty(opts.clusterColors)
            clusterColors = lines(K);
        else
            clusterColors = opts.clusterColors;
        end

        for k = 1:K
            range = groupAngleRanges(k,:); % [thetaStart thetaEnd]
            draw_cluster_wedge(range, opts.radius, clusterColors(k,:));
        end

        % Cluster name labels, placed at outer arc midpoint of each group
        if ~isempty(opts.clusterNames)
            for k = 1:K
                range = groupAngleRanges(k,:);
                midAngle = mean(range);
                lr = opts.radius * 1.45;
                text(lr*cos(midAngle), lr*sin(midAngle), opts.clusterNames{k}, ...
                    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                    'FontWeight','bold', 'FontSize', 12, 'Color', clusterColors(k,:)*0.7);
            end
        end
    end

    % --- Draw arrows for transitions (i ~= j) ---
    for i = 1:N
        for j = 1:N
            pij = P(i,j);
            if pij < opts.minProb
                continue;
            end

            if i == j
                if opts.selfLoop
                    draw_self_loop(x(i), y(i), theta_states(i), pij, opts);
                end
                continue;
            end

            draw_curved_arrow([x(i) y(i)], [x(j) y(j)], pij, opts);
        end
    end

    % --- Draw state nodes on top ---
    for i = 1:N
        plot(x(i), y(i), 'o', 'MarkerSize', sqrt(opts.nodeSize)*3, ...
            'MarkerFaceColor', [0.9 0.9 1], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        lx = x(i) * 1.18;
        ly = y(i) * 1.18;
        text(lx, ly, stateNames{i}, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize', 11);
    end

    if isempty(partition) || K == 1
        xlim([-1.5 1.5]*opts.radius);
        ylim([-1.5 1.5]*opts.radius);
    else
        xlim([-1.7 1.7]*opts.radius);
        ylim([-1.7 1.7]*opts.radius);
    end
    title('Transition Probability Diagram (Grouped by Cluster)');

    % --- Probability colorbar ---
    cb_ax = axes('Position',[0.92 0.3 0.02 0.4]);
    nshades = 256;
    cdata = zeros(nshades,1,3);
    for kk = 1:nshades
        cdata(kk,1,:) = opts.colormapFn((kk-1)/(nshades-1));
    end
    image(cb_ax, flipud(cdata));
    set(cb_ax,'XTick',[],'YTick',[0 nshades],'YTickLabel',{'1.0','0.0'}, ...
        'YAxisLocation','right','Box','on');
    ylabel(cb_ax,'P(i\rightarrowj)','Rotation',270,'VerticalAlignment','bottom');

end

% ===================== HELPER FUNCTIONS =====================

function theta = compute_even_angles(N)
    theta = linspace(0, 2*pi, N+1);
    theta(end) = [];
    theta = theta + pi/2;
    theta = -theta + pi;
end

function [theta_states, order, groupAngleRanges] = compute_grouped_angles(partIdx, gapFrac)
% Assigns angular positions to states such that states in the same
% cluster are contiguous, with gaps between clusters.
%
% theta_states     - N x 1 angles, in ORIGINAL state order (i.e.
%                     theta_states(i) is the angle for state i)
% order             - permutation giving the display order (cluster-major)
% groupAngleRanges  - K x 2, [thetaStart thetaEnd] for each cluster's span

    N = numel(partIdx);
    K = max(partIdx);

    % Order states cluster-major (stable within cluster)
    order = [];
    clusterSizes = zeros(K,1);
    for k = 1:K
        idx = find(partIdx == k);
        clusterSizes(k) = numel(idx);
        order = [order; idx]; %#ok<AGROW>
    end

    totalGap = 2*pi * gapFrac;
    if K > 1
        gapEach = totalGap / K;
    else
        gapEach = 0;
    end
    totalNodeAngle = 2*pi - K*gapEach;

    % Distribute node angle proportionally to cluster size
    theta_ordered = zeros(N,1);
    groupAngleRanges = zeros(K,2);

    angleCursor = pi/2; % start at top, go clockwise (matching original convention)
    for k = 1:K
        nK = clusterSizes(k);
        spanK = totalNodeAngle * (nK / N);

        startAngle = angleCursor;
        if nK == 1
            positions = startAngle - spanK/2; % center single node in its span
        else
            % Evenly space nK nodes across spanK, clockwise
            step = spanK / nK;
            positions = startAngle - (step*(0:nK-1) + step/2);
        end

        idxInOrder = sum(clusterSizes(1:k-1)) + (1:nK);
        theta_ordered(idxInOrder) = positions;

        groupAngleRanges(k,:) = sort([startAngle, startAngle - spanK]);

        angleCursor = startAngle - spanK - gapEach;
    end

    % Map back to original state indexing
    theta_states = zeros(N,1);
    theta_states(order) = theta_ordered;
end

function draw_cluster_wedge(angleRange, radius, color)
% Draws a faint filled wedge (annular sector) behind a cluster's nodes,
% as a visual grouping cue.

    thetaStart = angleRange(1);
    thetaEnd = angleRange(2);
    nArc = 20;
    thetas = linspace(thetaStart, thetaEnd, nArc);

    rOuter = radius * 1.32;
    rInner = radius * 0.0; % wedge from center

    xOuter = rOuter * cos(thetas);
    yOuter = rOuter * sin(thetas);

    xs = [0, xOuter, 0];
    ys = [0, yOuter, 0];

    fill(xs, ys, color, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end

function draw_curved_arrow(p1, p2, prob, opts)
    mid = (p1 + p2) / 2;
    dir = p2 - p1;
    perp = [-dir(2), dir(1)];
    perp = perp / norm(perp);

    curvature = 0.15 * norm(p1);
    ctrl = mid + perp * curvature;

    t = linspace(0,1,30)';
    bez = (1-t).^2 * p1 + 2*(1-t).*t * ctrl + t.^2 * p2;

    shrink = 0.12;
    bez = trim_curve_ends(bez, shrink);

    color = opts.colormapFn(prob);
    lw = 0.5 + 6*prob;

    plot(bez(:,1), bez(:,2), '-', 'Color', color, 'LineWidth', lw);

    tipDir = bez(end,:) - bez(end-1,:);
    tipDir = tipDir / norm(tipDir);
    arrow_len = 0.06;
    arrow_w = 0.03;
    tip = bez(end,:);
    base = tip - arrow_len * tipDir;
    perpTip = [-tipDir(2), tipDir(1)];
    p_left = base + arrow_w*perpTip;
    p_right = base - arrow_w*perpTip;

    fill([tip(1) p_left(1) p_right(1)], [tip(2) p_left(2) p_right(2)], ...
        color, 'EdgeColor', 'none');

    labelPt = bez(round(end/2),:);
    text(labelPt(1), labelPt(2), sprintf('%.2f', prob), ...
        'FontSize', 8, 'HorizontalAlignment','center', ...
        'BackgroundColor','w', 'Margin', 0.5);
end

function bez = trim_curve_ends(bez, shrinkFrac)
    n = size(bez,1);
    k = max(1, round(n*shrinkFrac));
    bez = bez(k+1:end-k, :);
end

function draw_self_loop(x, y, theta, prob, opts)
    loopR = 0.18;
    center = [x y] + loopR * [cos(theta), sin(theta)];

    ang = linspace(0.2, 2*pi-0.2, 30);
    loopPts = center + loopR*[cos(ang)' sin(ang)'];

    color = opts.colormapFn(prob);
    lw = 0.5 + 6*prob;

    plot(loopPts(:,1), loopPts(:,2), '-', 'Color', color, 'LineWidth', lw);

    tip = loopPts(end,:);
    tipDir = loopPts(end,:) - loopPts(end-1,:);
    tipDir = tipDir / norm(tipDir);
    arrow_len = 0.05;
    arrow_w = 0.025;
    base = tip - arrow_len*tipDir;
    perpTip = [-tipDir(2), tipDir(1)];
    p_left = base + arrow_w*perpTip;
    p_right = base - arrow_w*perpTip;
    fill([tip(1) p_left(1) p_right(1)], [tip(2) p_left(2) p_right(2)], ...
        color, 'EdgeColor', 'none');

    text(center(1), center(2), sprintf('%.2f', prob), ...
        'FontSize', 8, 'HorizontalAlignment','center', ...
        'BackgroundColor','w', 'Margin', 0.5);
end
