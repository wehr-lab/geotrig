function [glmResults] = pp_glm(eventFrames, eventnames, stateMask, statenames, ...
                                num_frames, fps, varargin)
% PP_GLM  Point process GLM: event rate ~ history of other events + current state
%   For each target event type, fits a Poisson GLM:
%     log(lambda_e(t)) = beta0 + sum_f sum_lag beta(f,lag)*history(f,t,lag)
%                               + sum_s gamma(s)*state(s,t)
%
%   [glmResults] = pp_glm(eventFrames, eventnames, stateMask, statenames, ...)
%
%   Optional name-value:
%     'HistoryLags' : lag bins in seconds (default [0.5 1 2 5])
%     'BinWidth'    : GLM bin width in seconds (default 0.5)
%     'Plot'        : true/false (default true)

p = inputParser;
p.addParameter('HistoryLags', [0.5 1 2 5]);
p.addParameter('BinWidth', 0.5);
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;

nE      = numel(eventnames);
nS      = numel(statenames);
binFr   = round(opt.BinWidth * fps);
nBins   = floor(num_frames / binFr);
lagsFr  = round(opt.HistoryLags * fps);
nLags   = numel(lagsFr);

% Build design matrix bins: [nBins x (nE*nLags + nS)]
% History features: for each event type f and lag l,
%   was there an event of type f in [t-lag_l, t-lag_{l-1}] ?
histFeats = zeros(nBins, nE * nLags);
for f = 1:nE
    ef = eventFrames{f};
    ef = ef(ef>=1 & ef<=num_frames);
    spike = zeros(num_frames,1); spike(ef) = 1;
    for l = 1:nLags
        lagEnd   = lagsFr(l);
        if l > 1
            lagStart = lagsFr(l-1);
        else
            lagStart = 0;
        end     
        col = (f-1)*nLags + l;
        for b = 1:nBins
            tCenter = (b-0.5)*binFr;
            idxLo = max(1, round(tCenter - lagEnd));
            idxHi = max(1, round(tCenter - lagStart));
            histFeats(b, col) = sum(spike(idxLo:idxHi));
        end
    end
end

% State features: fraction of bin spent in each state
stateFeat = zeros(nBins, nS);
for b = 1:nBins
    idx = ((b-1)*binFr+1) : min(b*binFr, num_frames);
    stateFeat(b,:) = mean(stateMask(idx,:), 1);
end

X = [histFeats, stateFeat];
featNames = {};
for f = 1:nE
    for l = 1:nLags
        featNames{end+1} = sprintf('%s lag%.1fs', eventnames{f}, opt.HistoryLags(l)); %#ok
    end
end
featNames = [featNames, statenames];

glmResults = struct();
cols = lines(nE);

if opt.Plot
    figure('Position',[100 100 1400 200*nE]);
    tiledlayout(nE, 1, 'TileSpacing','compact','Padding','compact');
end

for e = 1:nE
    % Response: spike count in each bin
    ef = eventFrames{e}; ef = ef(ef>=1 & ef<=num_frames);
    spike = zeros(num_frames,1); spike(ef)=1;
    y = zeros(nBins,1);
    for b=1:nBins
        idx = ((b-1)*binFr+1):min(b*binFr,num_frames);
        y(b) = sum(spike(idx));
    end

    % Fit Poisson GLM
    try
        mdl = fitglm(X, y, 'Distribution','poisson','Link','log', ...
                     'Intercept',true);
        coefs  = mdl.Coefficients.Estimate(2:end);
        pvals  = mdl.Coefficients.pValue(2:end);
        ci     = coefCI(mdl);
        ci     = ci(2:end,:);

        glmResults(e).event   = eventnames{e};
        glmResults(e).coefs   = coefs;
        glmResults(e).pvals   = pvals;
        glmResults(e).ci      = ci;
        glmResults(e).deviance = mdl.Deviance;
        glmResults(e).model   = mdl;

        if opt.Plot
            nexttile; hold on;
            nFeat = numel(coefs);
            sig = pvals < 0.05;
            errorbar(1:nFeat, coefs, coefs-ci(:,1), ci(:,2)-coefs, ...
                'o','Color',cols(e,:),'LineWidth',1.2,'MarkerSize',5);
            scatter(find(sig), coefs(sig), 60, 'r', 'filled');
            xline(nE*nLags+0.5,'k--');  % divider: history vs state features
            yline(0,'k:');
            xlim([0 nFeat+1]);
            set(gca,'XTick',1:nFeat,'XTickLabel',featNames,'XTickLabelRotation',40,...
                'FontSize',7,'TickLabelInterpreter','none');
            ylabel('log \beta (95% CI)');
            title(sprintf('GLM: %s  (dev=%.1f)', eventnames{e}, mdl.Deviance), ...
                'Interpreter','none','FontSize',9);
            text(nE*nLags*0.5, max(coefs)*0.9,'← history features', ...
                'FontSize',7,'HorizontalAlignment','center');
            text(nE*nLags + nS*0.5+0.5, max(coefs)*0.9,'state features →', ...
                'FontSize',7,'HorizontalAlignment','center');
        end

    catch ME
        warning('GLM failed for %s: %s', eventnames{e}, ME.message);
    end
end

if opt.Plot
    sgtitle('Point process GLM  [red dots = p<0.05, dashed line = history|state boundary]');
end
end