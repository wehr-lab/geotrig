function plot_unified_interaction(TPM, statenames, ...
                                   hawkesParams, eventnames, ...
                                   stateRates, ...
                                   glmResults, ...
                                   eventStateJS, ...
                                   varargin)
% PLOT_UNIFIED_INTERACTION  Unified state-event interaction diagram
%
% Inputs:
%   TPM          : [nS x nS] transition probability matrix
%   statenames   : cell [1 x nS]
%   hawkesParams : output of pp_hawkes_fit (uses .branchingRatio [nE x nE])
%   eventnames   : cell [1 x nE]
%   stateRates   : output of pp_state_rates (uses .foldChange [nS x nE])
%   glmResults   : output of pp_glm (uses .coefs, .pvals per event)
%   eventStateJS : [nE x nS] JS divergence of TPM row s after event e
%                  (from your event-conditioned TPM analysis, results(e).jsdiv)
%
% Optional name-value:
%   'TPM_thresh'       : min TPM value to draw (default 0.1)
%   'Hawkes_thresh'    : min branching ratio to draw (default 0.1)
%   'StateEvent_thresh': min log2 fold-change to draw state->event (default 1)
%   'EventState_thresh': min JS divergence to draw event->state (default 0.1)
%   'MaxLineWidth'     : max arrow line width (default 6)
%   'NodeRadius'       : node circle radius (default 0.07)
%   'RingRadiusState'  : radius of state ring (default 1.0)
%   'RingRadiusEvent'  : radius of event ring (default 0.52)

% Key design decisions:
% 
% Two rings cleanly separate the two node types while keeping cross-arrows
% short and readable. More legible than mixing everything on one ring. 
%
% Four arrow colors (gray / orange / blue / red) immediately tell you which
% of the four interaction types you're looking at — no need to read a
% complicated legend to parse each arrow.
%
% Thresholds are critical — start permissive and tighten until the diagram
% is readable. The *_thresh parameters are your main tool for decluttering.
%
% GLM gates the state→event arrows — fold-change alone can be misleading
% with small counts, so only arrows that are both large in fold-change and
% significant in the GLM are drawn.


% Usage:
% Assemble eventStateJS from your existing results struct:
% eventStateJS = zeros(nEvents, nStates);
% for e = 1:nEvents
%     eventStateJS(e,:) = results(e).jsdiv;
% end
% 
% plot_unified_interaction(TPM, statenames, hawkesParams, eventnames, ...
%     stateRates, glmResults, eventStateJS, ...
%     'TPM_thresh', 0.08, ...
%     'Hawkes_thresh', 0.08, ...
%     'StateEvent_thresh', 0.8, ...
%     'EventState_thresh', 0.08);


p = inputParser;
p.addParameter('TPM_thresh',        0.10);
p.addParameter('Hawkes_thresh',     0.10);
p.addParameter('StateEvent_thresh', 1.0);   % log2 fold-change
p.addParameter('EventState_thresh', 0.10);  % JS divergence
p.addParameter('MaxLineWidth',      6);
p.addParameter('NodeRadius',        0.07);
p.addParameter('RingRadiusState',   1.0);
p.addParameter('RingRadiusEvent',   0.52);
p.parse(varargin{:});
opt = p.Results;

nS = numel(statenames);
nE = numel(eventnames);

% ---- Node positions ----
% States: outer ring
thetaS = linspace(pi/2, pi/2 - 2*pi*(nS-1)/nS, nS);
xyS = opt.RingRadiusState * [cos(thetaS); sin(thetaS)]';

% Events: inner ring, rotated by half-step to interleave
thetaE = linspace(pi/2 + pi/nE, pi/2 + pi/nE - 2*pi*(nE-1)/nE, nE);
xyE = opt.RingRadiusEvent * [cos(thetaE); sin(thetaE)]';

% ---- Precompute max values for scaling ----
maxTPM    = max(TPM(:));
maxHawkes = max(hawkesParams.branchingRatio(:));
maxFC     = max(abs(log2(stateRates.foldChange(:) + eps)));
maxJS     = max(eventStateJS(:));

lw = @(val, maxVal) max(0.5, opt.MaxLineWidth * val / max(maxVal, eps));

figure('Position', [100 100 1000 1000]);
ax = axes; hold on; axis equal off;

%fontsize prefs
fs1=18; %labels
fs2=18;
fs3=24; %title

% ================================================================
% 1. STATE -> STATE arrows (TPM) — gray, outer ring self-referential
% ================================================================
for i = 1:nS
    for j = 1:nS
        if i==j, continue; end
        v = TPM(i,j);
        if v < opt.TPM_thresh, continue; end
        col = [1 1 1]*max(0, 0.7 - 0.7*v/maxTPM);  % darker = stronger
        draw_arrow(xyS(i,:), xyS(j,:), lw(v,maxTPM), col, opt.NodeRadius, 0.06, 0.10, ax);
    end
end

% ================================================================
% 2. EVENT -> EVENT arrows (Hawkes branching ratio) — orange, inner ring
% ================================================================
BR = hawkesParams.branchingRatio;
for e = 1:nE
    for f = 1:nE
        if e==f, continue; end
        v = BR(e,f);  % event f triggers event e
        if v < opt.Hawkes_thresh, continue; end
        col = [0.9 0.5 0.0] * min(1, 0.4 + 0.6*v/maxHawkes);
        draw_arrow(xyE(f,:), xyE(e,:), lw(v,maxHawkes), col, opt.NodeRadius*0.85, 0.05, 0.08, ax);
    end
end

% ================================================================
% 3. STATE -> EVENT arrows (rate enrichment) — blue, cross-ring
%    Use GLM significance as gate, fold-change as thickness
% ================================================================
FC = log2(stateRates.foldChange + eps);  % [nS x nE]
for s = 1:nS
    for e = 1:nE
        v = FC(s,e);
        % if (v) < opt.StateEvent_thresh, continue; end %mw
         if abs(v) < opt.StateEvent_thresh, continue; end
        % Check GLM significance if available
        if ~isempty(glmResults) && numel(glmResults) >= e
            % find the state feature index in GLM coefs
            nHistFeats = size(glmResults(e).coefs,1) - nS;
            stateCoefIdx = nHistFeats + s;
            if stateCoefIdx <= numel(glmResults(e).pvals)
                if glmResults(e).pvals(stateCoefIdx) > 0.05, continue; end
            end
        end
        if v > 0
            col = [0.1 0.3 0.9];  % blue = enriched
        else
            col = [0.4 0.7 1.0];  % light blue = depleted
        end
        draw_arrow_cross(xyS(s,:), xyE(e,:), lw(abs(v),maxFC), col, opt.NodeRadius, ax);
    end
end

% ================================================================
% 4. EVENT -> STATE arrows (post-event TPM perturbation) — red, cross-ring
%    eventStateJS: [nE x nS], JS divergence of row s after event e
% ================================================================
for e = 1:nE
    for s = 1:nS
        v = eventStateJS(e,s);
        if v < opt.EventState_thresh, continue; end
        intensity = min(1, 0.3 + 0.7*v/maxJS);
        col = [intensity 0.1 0.1];
        draw_arrow_cross(xyE(e,:), xyS(s,:), lw(v,maxJS), col, opt.NodeRadius, ax);
    end
end

% ================================================================
% Draw nodes (events first so states appear on top)
% ================================================================
% Event nodes — rounded squares, inner ring
for e = 1:nE
    r = opt.NodeRadius * 0.85;
    rectangle('Position', [xyE(e,1)-r, xyE(e,2)-r, 2*r, 2*r], ...
        'Curvature', 0.3, 'FaceColor', [1.0 0.92 0.75], ...
        'EdgeColor', [0.8 0.5 0.0], 'LineWidth', 1.5);
    labelPos = xyE(e,:) * (1 + 0.55/opt.RingRadiusEvent);
    text(labelPos(1), labelPos(2), eventnames{e}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize', fs1, 'Interpreter','none', 'FontWeight','bold', ...
        'Color', [0.6 0.3 0.0]);
end

% State nodes — circles, outer ring
for s = 1:nS
    r = opt.NodeRadius;
    rectangle('Position', [xyS(s,1)-r, xyS(s,2)-r, 2*r, 2*r], ...
        'Curvature', 1, 'FaceColor', [0.85 0.92 1.0], ...
        'EdgeColor', [0.1 0.2 0.7], 'LineWidth', 1.8);
    labelPos = xyS(s,:) * (1 + 0.45/opt.RingRadiusState);
    text(labelPos(1), labelPos(2), statenames{s}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize', fs2, 'Interpreter','none', 'FontWeight','bold', ...
        'Color', [0.05 0.15 0.55]);
end

% ================================================================
% Legend
% ================================================================
lx = -1.55; ly = -1.35; ldy = 0.10;
line([lx lx+0.12], [ly ly],     'Color',[0.5 0.5 0.5],'LineWidth',2); 
text(lx+0.15, ly,     'State\rightarrowState (TPM)',         'FontSize',fs1,'Interpreter','tex');
line([lx lx+0.12], [ly-ldy ly-ldy], 'Color',[0.9 0.5 0.0],'LineWidth',2);
text(lx+0.15, ly-ldy, 'Event\rightarrowEvent (Hawkes)',      'FontSize',fs1,'Interpreter','tex');
line([lx lx+0.12], [ly-2*ldy ly-2*ldy], 'Color',[0.1 0.3 0.9],'LineWidth',2);
text(lx+0.15, ly-2*ldy,'State\rightarrowEvent (rate enrich.)','FontSize',fs1,'Interpreter','tex');
line([lx lx+0.12], [ly-3*ldy ly-3*ldy], 'Color',[0.9 0.1 0.1],'LineWidth',2);
text(lx+0.15, ly-3*ldy,'Event\rightarrowState (TPM shift)',  'FontSize',fs1,'Interpreter','tex');

text(0, -1.6, 'Width = effect size   |   Threshold-filtered', ...
    'HorizontalAlignment','left','FontSize',fs1,'Color',[0.4 0.4 0.4]);

title('Unified State–Event Interaction Diagram', 'FontSize', fs3, 'FontWeight','bold');
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
end


% =========================================================================
% Helper: draw curved arrow between two nodes on the same ring
function draw_arrow(p1, p2, lw, col, nodeR, offsetAmt, arrowScale, ax)
    dir     = p2 - p1;
    L       = norm(dir);
    unitDir = dir / L;
    normal  = [-unitDir(2), unitDir(1)];

    startPt = p1 + unitDir*nodeR + normal*offsetAmt;
    endPt   = p2 - unitDir*nodeR + normal*offsetAmt;
    midPt   = (startPt+endPt)/2 + normal*0.12;

    t      = linspace(0,1,40);
    cx     = (1-t).^2*startPt(1) + 2*(1-t).*t*midPt(1) + t.^2*endPt(1);
    cy     = (1-t).^2*startPt(2) + 2*(1-t).*t*midPt(2) + t.^2*endPt(2);

    plot(ax, cx, cy, '-', 'Color', col, 'LineWidth', lw);
    draw_arrowhead(ax, cx, cy, col, lw, arrowScale);
end

% Helper: draw straight(ish) arrow between rings (cross connections)
function draw_arrow_cross(p1, p2, lw, col, nodeR, ax)
    dir     = p2 - p1;
    L       = norm(dir);
    unitDir = dir / L;

    startPt = p1 + unitDir * nodeR;
    endPt   = p2 - unitDir * nodeR;

    % Slight curve to reduce overlap with other cross arrows
    normal  = [-unitDir(2), unitDir(1)];
    midPt   = (startPt+endPt)/2 + normal*0.04;

    t  = linspace(0,1,30);
    cx = (1-t).^2*startPt(1) + 2*(1-t).*t*midPt(1) + t.^2*endPt(1);
    cy = (1-t).^2*startPt(2) + 2*(1-t).*t*midPt(2) + t.^2*endPt(2);

    plot(ax, cx, cy, '-', 'Color', col, 'LineWidth', lw);
    draw_arrowhead(ax, cx, cy, col, lw, 0.07);
end

% Helper: arrowhead scaled with line width
function draw_arrowhead(ax, cx, cy, col, lw, baseScale)
    arrowLen = baseScale * (0.6 + 0.08*lw);
    arrowWid = 0.45 + 0.04*lw;

    vec  = [cx(end)-cx(end-2), cy(end)-cy(end-2)];
    vec  = vec / norm(vec) * arrowLen;
    perp = [-vec(2), vec(1)] * arrowWid;
    base = [cx(end), cy(end)] - vec;

    fill(ax, [cx(end), base(1)+perp(1), base(1)-perp(1)], ...
             [cy(end), base(2)+perp(2), base(2)-perp(2)], ...
         col, 'EdgeColor','none');
end