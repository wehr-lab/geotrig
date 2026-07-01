function plot_tpm_diff_circle(diffMat, statenames, varargin)
% plot_tpm_diff_circle  Circular diagram of TPM differences
%   diffMat   : nStates x nStates matrix of (TPM_post - TPM_baseline)
%   statenames: cell array of state labels
%
% Optional name-value args:
%   'Title'      : title string
%   'MinAbsDiff' : threshold below which arrows are not drawn (default 0.01)
%   'MaxLineWidth': linewidth corresponding to max |diff| (default 8)
%   'NodeRadius' : radius of node circles (default 0.08)

p = inputParser;
p.addParameter('Title','');
p.addParameter('MinAbsDiff',0.01);
p.addParameter('MaxLineWidth',8);
p.addParameter('NodeRadius',0.08);
p.parse(varargin{:});
opt = p.Results;

nStates = numel(statenames);

% Node positions on unit circle
theta = linspace(pi/2, pi/2 - 2*pi*(nStates-1)/nStates, nStates);
xy = [cos(theta); sin(theta)]';

maxAbs = max(abs(diffMat(:)));
if maxAbs == 0, maxAbs = 1; end

hold on; axis equal off;

% --- Draw arrows for off-diagonal transitions ---
for i = 1:nStates
    for j = 1:nStates
        if i == j, continue; end
        d = diffMat(i,j);
        if abs(d) < opt.MinAbsDiff, continue; end

        p1 = xy(i,:);
        p2 = xy(j,:);

        % Shrink endpoints to node edge, offset curve slightly for i->j vs j->i
        dir = p2 - p1;
        L = norm(dir);
        unitDir = dir / L;
        normal = [-unitDir(2), unitDir(1)]; % perpendicular, for curvature/offset

        startPt = p1 + unitDir*opt.NodeRadius;
        endPt   = p2 - unitDir*opt.NodeRadius;

        % Offset so i->j and j->i don't overlap
        offsetAmt = 0.04;
        startPt = startPt + normal*offsetAmt;
        endPt   = endPt   + normal*offsetAmt;
        midPt   = (startPt+endPt)/2 + normal*0.08; % bow the curve outward

        % Bezier curve
        t = linspace(0,1,30);
        curveX = (1-t).^2*startPt(1) + 2*(1-t).*t*midPt(1) + t.^2*endPt(1);
        curveY = (1-t).^2*startPt(2) + 2*(1-t).*t*midPt(2) + t.^2*endPt(2);

        lw = max(0.5, opt.MaxLineWidth * abs(d)/maxAbs);
        if d > 0
            col = [0.85 0.1 0.1] * min(1, abs(d)/maxAbs + 0.3); % red, increase
            col = [1 0 0]*0.7 + [0 0 0]*0; %#ok<NASGU> 
            col = [1, 1-min(1,abs(d)/maxAbs), 1-min(1,abs(d)/maxAbs)]; % light->dark red
        else
            col = [1-min(1,abs(d)/maxAbs), 1-min(1,abs(d)/maxAbs), 1]; % light->dark blue
        end

        plot(curveX, curveY, '-', 'Color', col, 'LineWidth', lw);

        % % Arrowhead at end
        % arrowVec = [curveX(end)-curveX(end-1), curveY(end)-curveY(end-1)];
        % arrowVec = arrowVec / norm(arrowVec) * 0.04;
        % perp = [-arrowVec(2), arrowVec(1)]*0.6;
        % ah = [curveX(end), curveY(end)] - arrowVec;
        % fill([curveX(end), ah(1)+perp(1), ah(1)-perp(1)], ...
        %     [curveY(end), ah(2)+perp(2), ah(2)-perp(2)], col, 'EdgeColor','none');

        % Arrowhead at end — scale with line width
        arrowLen = 0.025 + 0.012*lw;   % length of arrowhead, grows with line weight
        arrowWidth = 0.4 + 0.05*lw;    % half-width factor, grows with line weight

        arrowVec = [curveX(end)-curveX(end-1), curveY(end)-curveY(end-1)];
        arrowVec = arrowVec / norm(arrowVec) * arrowLen;
        perp = [-arrowVec(2), arrowVec(1)] * arrowWidth;
        ah = [curveX(end), curveY(end)] - arrowVec;
        fill([curveX(end), ah(1)+perp(1), ah(1)-perp(1)], ...
            [curveY(end), ah(2)+perp(2), ah(2)-perp(2)], col, 'EdgeColor','none');
    end
end

% --- Self-loops on diagonal ---
for i = 1:nStates
    d = diffMat(i,i);
    if abs(d) < opt.MinAbsDiff, continue; end
    center = xy(i,:);
    outward = center / norm(center); % direction away from circle center
    loopCenter = center + outward*0.18;
    th = linspace(0, 2*pi, 30);
    r = 0.06;
    lx = loopCenter(1) + r*cos(th);
    ly = loopCenter(2) + r*sin(th);

    lw = max(0.5, opt.MaxLineWidth * abs(d)/maxAbs);
    if d > 0
        col = [1, 1-min(1,abs(d)/maxAbs), 1-min(1,abs(d)/maxAbs)];
    else
        col = [1-min(1,abs(d)/maxAbs), 1-min(1,abs(d)/maxAbs), 1];
    end
    plot(lx, ly, '-', 'Color', col, 'LineWidth', lw);
end

% --- Draw nodes ---
for i = 1:nStates
    rectangle('Position',[xy(i,1)-opt.NodeRadius, xy(i,2)-opt.NodeRadius, ...
               2*opt.NodeRadius, 2*opt.NodeRadius], ...
               'Curvature',1, 'FaceColor','w', 'EdgeColor','k', 'LineWidth',1.2);
    labelPos = xy(i,:) * 1.28;
    text(labelPos(1), labelPos(2), statenames{i}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',10, 'interpreter', 'none');
end

xlim([-1.5 1.5]); ylim([-1.5 1.5]);
title(opt.Title, 'Interpreter','none');

% Legend / colorbar proxy
text(1.3, -1.4, 'red = \uparrow prob, blue = \downarrow prob, width = |\Deltap|', ...
    'FontSize',8, 'HorizontalAlignment','right');
end