function [stateRates, rateMat] = pp_state_rates(eventFrames, eventnames, ...
                                                 stateMask, statenames, fps, varargin)
% PP_STATE_RATES  Event rates conditioned on current behavioral state
%   Answers: is event X more likely during state Y?
%
%   [stateRates, rateMat] = pp_state_rates(...)
%
%   Optional name-value:
%     'Plot' : true/false (default true)

p = inputParser;
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;

nEvents = numel(eventnames);
nStates = size(stateMask, 2);
num_frames = size(stateMask, 1);

rateMat = zeros(nStates, nEvents);  % [state x event] in events/s

for s = 1:nStates
    stateDur = sum(stateMask(:,s)) / fps;  % total time in state (s)
    if stateDur == 0, continue; end
    for e = 1:nEvents
        ef = eventFrames{e};
        ef = ef(ef >= 1 & ef <= num_frames);
        % Count events that occurred while in this state
        nOccurred = sum(stateMask(ef, s));
        rateMat(s,e) = nOccurred / stateDur;
    end
end

% Also compute overall (marginal) rate for each event
overallRate = zeros(1, nEvents);
for e = 1:nEvents
    ef = eventFrames{e};
    overallRate(e) = numel(ef) / (num_frames/fps);
end

% Fold-change relative to overall rate
foldChange = rateMat ./ (overallRate + eps);

stateRates.rateMat    = rateMat;
stateRates.foldChange = foldChange;
stateRates.overallRate = overallRate;

if opt.Plot
    figure('Position',[100 100 1200 430]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(rateMat);
    colorbar; colormap(hot);
    set(gca,'XTick',1:nEvents,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nStates,'YTickLabel',statenames,'TickLabelInterpreter','none');
    title('Event rate (events/s) by state');
    for s=1:nStates
        for e=1:nEvents
            text(e,s,sprintf('%.3f',rateMat(s,e)),'HorizontalAlignment','center',...
                'FontSize',7,'Color','w');
        end
    end

    nexttile;
    cmax = max(abs(log2(foldChange(:)+eps)));
    imagesc(log2(foldChange+eps), [-cmax cmax]);
    colorbar; colormap(gca, redblue_cmap());
    set(gca,'XTick',1:nEvents,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nStates,'YTickLabel',statenames,'TickLabelInterpreter','none');
    title('log_2 fold-change vs. overall rate  (red=enriched, blue=depleted)');
    for s=1:nStates
        for e=1:nEvents
            text(e,s,sprintf('%.1fx',foldChange(s,e)),'HorizontalAlignment','center',...
                'FontSize',7);
        end
    end
end
end