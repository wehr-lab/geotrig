function plot_state_trajectories(TPM, statenames, varargin)
% plot_state_trajectories  Plot P(state at t | start state) vs time
%   TPM        : [nStates x nStates] row-stochastic
%   statenames : cell array of state name strings
%
% Optional name-value args:
%   'MaxSteps' : number of steps to simulate (default 100)
%   'fps'      : frames per second for time axis (default 200)
%   'StartStates' : indices of starting states to plot (default all)

p = inputParser;
p.addParameter('MaxSteps', 100);
p.addParameter('fps', 200);
p.addParameter('StartStates', 1:size(TPM,1));
p.parse(varargin{:});
opt = p.Results;

nStates  = size(TPM, 1);
maxSteps = opt.MaxSteps;
tAxis    = (1:maxSteps) / opt.fps;
tLabel   = opt.fps > 1;  % show seconds if fps provided

% Precompute TPM^t for all t
TPMpow = zeros(nStates, nStates, maxSteps);
TPMpow(:,:,1) = TPM;
for t = 2:maxSteps
    TPMpow(:,:,t) = TPMpow(:,:,t-1) * TPM;
end

cols      = parula(nStates);
startList = opt.StartStates;
nPlots    = numel(startList);
nCols     = 2;
nRows     = ceil(nPlots / nCols);

figure('Position', [100 100 1200 220*nRows]);
tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

for k = 1:nPlots
    i = startList(k);
    nexttile; hold on;
    for j = 1:nStates
        plot(tAxis, squeeze(TPMpow(i,j,:)), '-', 'Color', cols(j,:), 'LineWidth', 1.4);
    end
    title(sprintf('Start: %s', statenames{i}), 'Interpreter','none', 'FontSize', 9);
    if tLabel
        xlabel('Time (s)');
    else
        xlabel('Steps');
    end
    ylabel('Probability');
    ylim([0 1]); xlim([tAxis(1) tAxis(end)]);
    grid on;
end

% Single shared legend on last tile
nexttile(nRows*nCols);  % may reuse last tile if nPlots is even — use a separate axes instead
lgd = legend(statenames, 'Interpreter','none', 'FontSize', 8, 'Location', 'best');

sgtitle('P(state at time t | starting state)', 'FontSize', 12);
end