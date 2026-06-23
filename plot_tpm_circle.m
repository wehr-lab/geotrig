function plot_tpm_circle(P, stateNames, opts)
% PLOT_TPM_CIRCLE  Visualize a transition probability matrix (TPM) as a
% circular diagram, with arrows whose thickness/darkness encode the
% transition probability P(i,j).
%
% USAGE:
%   plot_tpm_circle(P)
%   plot_tpm_circle(P, stateNames)
%   plot_tpm_circle(P, stateNames, opts)
%
% INPUTS:
%   P          - N x N transition probability matrix (rows sum to 1)
%   stateNames - (optional) cell array of N strings, e.g. {'A','B','C'}
%   opts       - (optional) struct with optional fields:
%       .minProb     - probabilities below this are not drawn (default 0.01)
%       .radius      - circle radius for state placement (default 1)
%       .nodeSize    - marker size for state circles (default 60)
%       .colormapFn  - function handle for arrow color, e.g. @gray (default)
%       .selfLoop    - true/false, draw self-transitions as loops (default true)
%       .title       - title (default 'Transition Probability Diagram')
%       .colorbar    - show a colorbar for probability scale (default 1)
%
% EXAMPLE:
%   P = [0.7 0.2 0.1;
%        0.3 0.4 0.3;
%        0.1 0.1 0.8];
%   plot_tpm_circle(P, {'Sunny','Cloudy','Rainy'});

    N = size(P,1);

    if nargin < 2 || isempty(stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    end
    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts,'minProb'),   opts.minProb   = 0.01; end
    if ~isfield(opts,'radius'),    opts.radius    = 1;    end
    if ~isfield(opts,'nodeSize'),  opts.nodeSize  = 60;   end
    if ~isfield(opts,'colormapFn'),opts.colormapFn= @(x) (1-x)*[1 1 1]; end % white->black
    if ~isfield(opts,'selfLoop'),  opts.selfLoop  = true; end
    if ~isfield(opts,'title'),  opts.title  = 'Transition Probability Diagram'; end
    if ~isfield(opts,'colorbar'),  opts.colorbar  = 1; end
    
   
    % --- Compute node positions evenly around a circle ---
    theta = linspace(0, 2*pi, N+1);
    theta(end) = [];
    theta = theta + pi/2;          % start at top
    theta = -theta + pi;           % go clockwise
    x = opts.radius * cos(theta);
    y = opts.radius * sin(theta);

    figure('Color','w');
    hold on; axis equal off;

    % --- Draw arrows for transitions (i ~= j) ---
    for i = 1:N
        for j = 1:N
            p = P(i,j);
            if p < opts.minProb
                continue;
            end

            if i == j
                if opts.selfLoop
                    draw_self_loop(x(i), y(i), theta(i), p, opts);
                end
                continue;
            end

            draw_curved_arrow([x(i) y(i)], [x(j) y(j)], p, opts);
        end
    end

    % --- Draw state nodes on top ---
    for i = 1:N
        plot(x(i), y(i), 'o', 'MarkerSize', sqrt(opts.nodeSize)*3, ...
            'MarkerFaceColor', [0.9 0.9 1], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        % Place label slightly outside the circle so it doesn't overlap arrows
        lx = x(i) * 1.18;
        ly = y(i) * 1.18;
        text(lx, ly, stateNames{i}, 'HorizontalAlignment','center', 'interpreter', 'none',...
            'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize', 11);
    end

    xlim([-1.5 1.5]*opts.radius);
    ylim([-1.5 1.5]*opts.radius);
    title(opts.title);

    % --- Add a simple legend/colorbar for probability scale ---
    if opts.colorbar
        cb_ax = axes('Position',[0.92 0.3 0.02 0.4]);
        nshades = 256;
        cdata = zeros(nshades,1,3);
        for k = 1:nshades
            cdata(k,1,:) = opts.colormapFn((k-1)/(nshades-1));
        end
        image(cb_ax, flipud(cdata));
        set(cb_ax,'XTick',[],'YTick',[0 nshades],'YTickLabel',{'1.0','0.0'}, ...
            'YAxisLocation','right','Box','on');
        ylabel(cb_ax,'P(i\rightarrowj)','Rotation',270,'VerticalAlignment','bottom');
    end
end

% ===================== HELPER FUNCTIONS =====================

function draw_curved_arrow(p1, p2, prob, opts)
% Draws a curved arrow from p1 to p2, bowing outward slightly so that
% bidirectional transitions (i->j and j->i) don't overlap.

    % Midpoint and perpendicular offset for curvature
    mid = (p1 + p2) / 2;
    dir = p2 - p1;
    perp = [-dir(2), dir(1)];
    perp = perp / norm(perp);

    curvature = 0.15 * norm(p1); % scale offset relative to circle size
    ctrl = mid + perp * curvature;

    % Quadratic Bezier curve points
    t = linspace(0,1,30)';
    bez = (1-t).^2 * p1 + 2*(1-t).*t * ctrl + t.^2 * p2;

    % Shorten the curve slightly at both ends so arrow doesn't overlap node markers
    shrink = 0.12;
    bez = trim_curve_ends(bez, shrink);

    color = opts.colormapFn(prob);
    lw = 0.5 + 6*prob; % line width scales with probability

    plot(bez(:,1), bez(:,2), '-', 'Color', color, 'LineWidth', lw);

    % Draw arrowhead at the end
    tipDir = bez(end,:) - bez(end-1,:);
    tipDir = tipDir / norm(tipDir);
    % arrow_len = 0.06;
    % arrow_w = 0.03;
    % Arrowhead at end — scale with line width
    arrow_len = 0*0.025 + 0.03*lw;   % length of arrowhead, grows with line weight
    arrow_w = 0*0.4 + 0.02*lw;    % half-width factor, grows with line weight

    tip = bez(end,:);
    base = tip - arrow_len * tipDir;
    perpTip = [-tipDir(2), tipDir(1)];
    p_left = base + arrow_w*perpTip;
    p_right = base - arrow_w*perpTip;

    fill([tip(1) p_left(1) p_right(1)], [tip(2) p_left(2) p_right(2)], ...
        color, 'EdgeColor', 'none');

    % Probability label near midpoint of curve
    labelPt = bez(round(end/2),:);
    text(labelPt(1), labelPt(2), sprintf('%.2f', prob), ...
        'FontSize', 8, 'HorizontalAlignment','center', ...
        'BackgroundColor','w', 'Margin', 0.5);
end

function bez = trim_curve_ends(bez, shrinkFrac)
% Removes points near both ends of the curve (fraction of total points)
    n = size(bez,1);
    k = max(1, round(n*shrinkFrac));
    bez = bez(k+1:end-k, :);
end

function draw_self_loop(x, y, theta, prob, opts)
% Draws a small loop above/outside the node to represent P(i,i)

    loopR = 0.18;
    center = [x y] + loopR * [cos(theta), sin(theta)];

    ang = linspace(0.2, 2*pi-0.2, 30);
    loopPts = center + loopR*[cos(ang)' sin(ang)'];

    color = opts.colormapFn(prob);
    lw = 0.5 + 6*prob;

    plot(loopPts(:,1), loopPts(:,2), '-', 'Color', color, 'LineWidth', lw);

    % Arrowhead at end of loop
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