% =========================================================================
%  POINT PROCESS ANALYSIS SUITE
%  Behavioral event & state analysis
%
%  Inputs (expected in workspace or passed to each function):
%    statenames  : {'hot pursuit','chase','following','stalk','wander','pause'}
%    eventnames  : {'failed_approach','target_loss','contact_loss',
%                   'contact_gain','intercept','cricket_jump','rangemin'}
%    eventFrames : cell array [1 x nEvents], each cell = column vector of
%                  frame indices on which that event occurred
%    stateMask   : [num_frames x nStates] logical, mutually exclusive
%    fps         : frames per second
%
% =========================================================================

% /Applications/MATLAB_R2023b.app/bin/matlab  -nodisplay -nodesktop -sd /Users/wehr/Documents/Analysis/geotrig -batch "pp_analysis_suite" -logfile /Users/wehr/Documents/Analysis/geotrig/pplog.txt

% save /Users/wehr/Documents/Analysis/geotrig/datacache num_frames   failed_approach_event_frames  target_loss_event_frames contact_loss_event_frames  ...
%     contact_gain_event_frames intercept_event_frames cricket_jump_event_frames rangemin_event_frames  ...
%     hotpursuit chase follow stalk  wander  pause;

if exist('hotpursuit')~=1
    load ~/Documents/Analysis/geotrig/datacache
end

close all

statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
eventnames = {'failed_approach','contact_loss', ...
    'contact_gain','intercept','cricket_jump','rangemin'};

nStates = numel(statenames);
nEvents = numel(eventnames);

eventFrames={failed_approach_event_frames,contact_loss_event_frames, ...
    contact_gain_event_frames,intercept_event_frames,cricket_jump_event_frames,rangemin_event_frames};
stateMask = [    hotpursuit(:),chase(:),follow(:),...
    stalk(:), wander(:), pause(:)];
fps=200;

pp_summary(eventFrames, eventnames, stateMask, statenames, num_frames, fps)
return

pp_state_rates(eventFrames, eventnames, stateMask, statenames, fps);


pp_crossintensity(eventFrames, eventnames, num_frames, fps);

%run hawkes fit and similate, which are not run by pp_summary
  hawkes_params = pp_hawkes_fit(eventFrames, eventnames, num_frames, fps);
  save /Users/wehr/Documents/Analysis/geotrig/datacache hawkes_params -append
% duration=100; %seconds
% simTimes = pp_hawkes_simulate(hawkes_params, eventnames, duration);
pp_hawkes_plot( eventnames, hawkes_params)


%print to pdf
f=findobj('type', 'figure');
pdffilename=sprintf('~/Documents/Analysis/geotrig/pp_analysis-figs-%s.pdf', replace(datestr(now), whitespacePattern, '-'));
fprintf('\nprinting %d figures to %s\n', length(f), pdffilename)
for idx=1:length(f)
    builtin('pause',.5)
    exportgraphics(f(idx),pdffilename,'Append',true)
    builtin('pause',.5)
    fprintf('.')
end
fprintf('done\n')

% =========================================================================
%  FUNCTION INDEX
%  1.  pp_rate             — kernel-smoothed event rate over session
%  2.  pp_iei              — inter-event interval distributions
%  3.  pp_crossintensity   — cross-intensity functions between event pairs
%  4.  pp_state_rates      — event rates conditioned on behavioral state
%  5.  pp_hawkes_fit       — fit a Hawkes process model (MLE)
%  6.  pp_hawkes_simulate  — simulate from fitted Hawkes model
%  7.  pp_glm              — point process GLM (event rate ~ history + state)
%  8.  pp_summary          — run and plot all of the above
% =========================================================================


%% =========================================================================
function [rate, tAxis] = pp_rate(eventFrames, eventnames, num_frames, fps, varargin)
% PP_RATE  Kernel-smoothed instantaneous event rate for each event type
%
%   [rate, tAxis] = pp_rate(eventFrames, eventnames, num_frames, fps)
%
%   Optional name-value:
%     'KernelWidth'  : smoothing kernel SD in seconds (default 1s)
%     'KernelType'   : 'gaussian' (default) or 'causal' (one-sided exponential)
%     'Plot'         : true/false (default true)

p = inputParser;
p.addParameter('KernelWidth', 1.0);   % seconds
p.addParameter('KernelType', 'gaussian');
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;

nEvents = numel(eventnames);
tAxis   = (1:num_frames)' / fps;
rate    = zeros(num_frames, nEvents);

kwFrames = round(opt.KernelWidth * fps);
halfW    = 3 * kwFrames;
t_kern   = (-halfW:halfW)';

for e = 1:nEvents
    % Binary spike train
    spike = zeros(num_frames, 1);
    ef = eventFrames{e};
    ef = ef(ef >= 1 & ef <= num_frames);
    spike(ef) = 1;

    % Build kernel
    switch opt.KernelType
        case 'gaussian'
            kern = exp(-0.5*(t_kern/kwFrames).^2);
        case 'causal'
            kern = zeros(size(t_kern));
            kern(t_kern >= 0) = exp(-t_kern(t_kern >= 0)/kwFrames);
    end
    kern = kern / sum(kern) * fps;  % normalize to units of events/s

    rate(:,e) = conv(spike, kern, 'same');
end

if opt.Plot
    cols = lines(nEvents);
    figure('Position',[100 100 1400 600]);
    tiledlayout(nEvents, 1, 'TileSpacing','compact','Padding','compact');
    for e = 1:nEvents
        nexttile; hold on;
        plot(tAxis, rate(:,e), '-', 'Color', cols(e,:), 'LineWidth', 1.2);
        % Rug plot of actual events
        ef = eventFrames{e} / fps;
        plot(ef, zeros(size(ef)), '|', 'Color', cols(e,:)*0.6, 'MarkerSize', 6);
        ylabel('Rate (Hz)', 'FontSize', 7);
        title(eventnames{e}, 'Interpreter','none', 'FontSize', 8);
        xlim([0 num_frames/fps]);
        if e < nEvents, set(gca,'XTickLabel',[]); end
    end
    xlabel('Time (s)');
    sgtitle(sprintf('Kernel-smoothed event rates  (kernel: %s, \\sigma=%.1fs)', ...
        opt.KernelType, opt.KernelWidth));
end
end


%% =========================================================================
function pp_iei(eventFrames, eventnames, fps, varargin)
% PP_IEI  Inter-event interval distributions for each event type
%         Tests for Poisson (exponential IEI) vs. clustered/refractory
%
%   pp_iei(eventFrames, eventnames, fps)
%
%   Optional name-value:
%     'MaxIEI'    : max IEI to display in seconds (default 30)
%     'nBins'     : histogram bins (default 40)

p = inputParser;
p.addParameter('MaxIEI', 30);
p.addParameter('nBins', 40);
p.parse(varargin{:});
opt = p.Results;

nEvents = numel(eventnames);
cols    = lines(nEvents);

figure('Position',[100 100 1400 700]);
tiledlayout(ceil(nEvents/2), 4, 'TileSpacing','compact','Padding','compact');

for e = 1:nEvents
    ef  = sort(eventFrames{e}(:)) / fps;
    iei = diff(ef);
    iei = iei(iei <= opt.MaxIEI & iei > 0);

    if numel(iei) < 3
        continue;
    end

    % Histogram
    nexttile; hold on;
    h = histogram(iei, opt.nBins, 'FaceColor', cols(e,:), ...
        'EdgeColor','none','Normalization','pdf');
    % Fit exponential (Poisson null)
    muHat = mean(iei);
    xFit  = linspace(0, opt.MaxIEI, 200);
    plot(xFit, exppdf(xFit, muHat), 'k--', 'LineWidth', 1.5);
    xlabel('IEI (s)'); ylabel('Density');
    title(eventnames{e}, 'Interpreter','none','FontSize',8);
    legend('data',sprintf('Exp(\\mu=%.1fs)',muHat),'Location','northeast','FontSize',6);

    % Coefficient of variation (CV=1 → Poisson, <1 → regular, >1 → bursty)
    cv = std(iei)/mean(iei);
    text(0.98, 0.95, sprintf('CV=%.2f\nn=%d', cv, numel(iei)), ...
        'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','top','FontSize',7);

    % Survival function (log scale reveals exponential = Poisson)
    nexttile; hold on;
    [f, x] = ecdf(iei);
    semilogy(x, 1-f, '-', 'Color', cols(e,:), 'LineWidth', 1.5);
    semilogy(xFit, expcdf(xFit, muHat, 'upper'), 'k--', 'LineWidth', 1.5);
    xlabel('IEI (s)'); ylabel('P(IEI > t)');
    title(sprintf('%s — survival', eventnames{e}), 'Interpreter','none','FontSize',8);
    grid on;
    ylim([1e-3 1]);
end
sgtitle('Inter-event interval distributions  (dashed = Poisson fit, CV: 1=Poisson, <1=regular, >1=bursty)');
end


%% =========================================================================
function [xcf, lags] = pp_crossintensity(eventFrames, eventnames, num_frames, fps, varargin)
% PP_CROSSINTENSITY  Pairwise cross-intensity (cross-correlogram) between events
%   Answers: does event A increase/decrease rate of event B in following seconds?
%
%   [xcf, lags] = pp_crossintensity(eventFrames, eventnames, num_frames, fps)
%
%   Optional name-value:
%     'MaxLag'       : max lag in seconds (default 5)
%     'BinWidth'     : bin width in seconds (default 0.2)
%     'nShuffles'    : number of shuffles for null distribution (default 200)
%     'Plot'         : true/false (default true)

p = inputParser;
p.addParameter('MaxLag', 5);
p.addParameter('BinWidth', 0.2);
p.addParameter('nShuffles', 200);
% p.addParameter('nShuffles', 3);
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;

nEvents  = numel(eventnames);
binW     = round(opt.BinWidth * fps);
maxLagFr = round(opt.MaxLag * fps);
lagBins  = -maxLagFr:binW:maxLagFr;
lags     = lagBins(1:end-1)/fps + opt.BinWidth/2;
nLags    = numel(lags);

xcf      = zeros(nEvents, nEvents, nLags);
xcf_null = zeros(nEvents, nEvents, nLags, opt.nShuffles);

% Convert to binary frame vectors
spikes = false(num_frames, nEvents);
for e = 1:nEvents
    ef = eventFrames{e};
    ef = ef(ef >= 1 & ef <= num_frames);
    spikes(ef, e) = true;
end

fprintf('\n computing cross-intensity...');
nbytes=  fprintf('%d/%d', 0, nEvents);
for ref = 1:nEvents
    fprintf(repmat('\b',1,nbytes));
    nbytes=  fprintf('%d/%d', ref, nEvents);
    refTimes = find(spikes(:,ref));
    if isempty(refTimes), continue; end

    for tgt = 1:nEvents
        tgtSpikes = spikes(:,tgt);

        % Accumulate peri-event histogram
        peh = zeros(nLags, 1);
        for k = 1:numel(lagBins)-1
            lo = lagBins(k); hi = lagBins(k+1);
            cnt = 0;
            for rt = refTimes'
                idx = (rt+lo):(rt+hi-1);
                idx = idx(idx>=1 & idx<=num_frames);
                cnt = cnt + sum(tgtSpikes(idx));
            end
            peh(k) = cnt;
        end
        % Normalize to rate (events/s)
        xcf(ref,tgt,:) = peh / (numel(refTimes) * opt.BinWidth);

        % Shuffle null: circularly shift target spike train
        shifts = randi(num_frames, opt.nShuffles, 1);
        for s = 1:opt.nShuffles
            tgtShuf = circshift(tgtSpikes, shifts(s));
            peh_s = zeros(nLags,1);
            for k = 1:numel(lagBins)-1
                lo = lagBins(k); hi = lagBins(k+1);
                cnt = 0;
                for rt = refTimes'
                    idx = (rt+lo):(rt+hi-1);
                    idx = idx(idx>=1 & idx<=num_frames);
                    cnt = cnt + sum(tgtShuf(idx));
                end
                peh_s(k) = cnt;
            end
            xcf_null(ref,tgt,:,s) = peh_s / (numel(refTimes) * opt.BinWidth);
        end
    end
end

if opt.Plot
    null_lo = prctile(xcf_null, 2.5, 4);
    null_hi = prctile(xcf_null, 97.5, 4);

    cols = lines(nEvents);
    figure('Position',[100 100 300*nEvents 220*nEvents]);
    tiledlayout(nEvents, nEvents, 'TileSpacing','compact','Padding','none');
    for ref = 1:nEvents
        for tgt = 1:nEvents
            nexttile; hold on;
            fill([lags, fliplr(lags)], ...
                 [squeeze(null_lo(ref,tgt,:))', fliplr(squeeze(null_hi(ref,tgt,:))')], ...
                 [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.6);
            plot(lags, squeeze(xcf(ref,tgt,:)), '-', 'Color', cols(tgt,:), 'LineWidth',1.2);
            xline(0, 'k--','LineWidth',0.8);
            if ref==tgt, set(gca,'Color',[0.97 0.97 0.93]); end  % highlight auto
            xlim([-opt.MaxLag opt.MaxLag]);
            set(gca,'FontSize',18,'XTick',[-opt.MaxLag 0 opt.MaxLag]);
            if tgt==1, ylabel(eventnames{ref},'Interpreter','none','FontSize',18,'Rotation',30,'HorizontalAlignment','right'); end
            if ref==1, title(eventnames{tgt},'Interpreter','none','FontSize',18); end
        end
    end
    sgtitle('Cross-intensity functions  [row=reference, col=target, lag=t(col)-t(row)]  gray=95% shuffle CI');
end
end


%% =========================================================================
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


%% =========================================================================
function  pp_hawkes_plot( eventnames, params)
%just the plot from hawkes_fit without the fitting
nE = numel(eventnames);
    figure('Position',[100 100 1300 430]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(params.alpha); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('\alpha  (excitation amplitude)'); colormap(hot);

    nexttile;
    imagesc(1./params.beta); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('1/\beta  (decay timescale, s)'); colormap(hot);

    nexttile;
    br = params.branchingRatio;
    imagesc(br, [0 1]); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('\alpha/\beta  (branching ratio: fraction of triggered events)');
    colormap(gca, hot);
    xlabel('triggering event')
    ylabel('received input')
    

    rho=max(abs(eig(params.branchingRatio)));
    sgtitle(sprintf('Hawkes process fit  [row=triggered event, col=triggering event], rho =%.2f', rho));
fprintf('\ncriticality: rho =%.4f', rho)
fprintf(['\nρ(K) < 1 → the process is stationary; cascades die out', ...
'\nρ(K) = 1 → critical; cascades can persist indefinitely (unit root behavior)', ...
'\nρ(K) > 1 → explosive; the process is non-stationary']);

opts.minProb=.2;
opts.title='Event interactions (Hawkes model)';
plot_tpm_circle(params.branchingRatio', eventnames, opts)

end

%% =========================================================================
function [params] = pp_hawkes_fit(eventFrames, eventnames, num_frames, fps, varargin)
% PP_HAWKES_FIT  Fit a multivariate Hawkes process via MLE
%   Each event type has:
%     mu(e)         : baseline rate (events/s)
%     alpha(e,f)    : excitation amplitude from event f onto event e
%     beta(e,f)     : decay rate of excitation from f onto e
%
%   Model: lambda_e(t) = mu(e) + sum_f sum_{t_f < t} alpha(e,f)*exp(-beta(e,f)*(t-t_f))
%
%   [params] = pp_hawkes_fit(eventFrames, eventnames, num_frames, fps)
%
%   Optional name-value:
%     'MaxDecay'  : max allowed decay time constant in seconds (default 10)
%     'Plot'      : true/false (default true)

p = inputParser;
p.addParameter('MaxDecay', 10);
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;  

fprintf('\n starting hawkes fit at %s', datestr(now))
htic=tic;

nE = numel(eventnames);
T  = num_frames / fps;  % total duration in seconds

% Convert eventFrames to times in seconds
eventTimes = cellfun(@(ef) sort(ef(:))/fps, eventFrames, 'UniformOutput',false);

% Parameter vector layout:
%   [mu(1..nE), alpha(1..nE*nE), beta(1..nE*nE)]
%   total: nE + 2*nE^2 parameters
nParams = nE + 2*nE^2;

% Initial guess: small baseline, weak excitation, moderate decay
mu0    = cellfun(@numel, eventTimes) / T;
alpha0 = 0.1 * ones(nE);
beta0  = ones(nE);   % 1 s decay

x0 = [mu0, alpha0(:)', beta0(:)'];

% Bounds: all positive, alpha < 1 (stationarity), beta > 0.1
lb = [zeros(1,nE),   zeros(1,nE^2),  0.1*ones(1,nE^2)];
ub = [10*ones(1,nE), ones(1,nE^2),   (1/opt.MaxDecay)*ones(1,nE^2)*10];

opts = optimoptions('fmincon','Display','off','MaxIterations',500, ...
                    'OptimalityTolerance',1e-6,'Algorithm','interior-point');
% opts = optimoptions('fmincon','Display','iter','MaxIterations',10, ...
%                     'OptimalityTolerance',1e-4,'Algorithm','interior-point');

negLogLik = @(x) hawkes_nll(x, eventTimes, nE, T);
[xOpt, fval] = fmincon(negLogLik, x0, [], [], [], [], lb, ub, [], opts);

% Unpack
params.mu    = xOpt(1:nE);
params.alpha = reshape(xOpt(nE+1:nE+nE^2),          nE, nE);
params.beta  = reshape(xOpt(nE+nE^2+1:nE+2*nE^2),   nE, nE);
params.nll   = fval;
params.BIC   = 2*fval + nParams*log(sum(cellfun(@numel,eventTimes)));
params.branchingRatio = params.alpha ./ params.beta;  % fraction of events that are "children"

fprintf('\n finished hawkes fit at %s', datestr(now))
d = seconds(toc(htic));
d.Format = 'hh:mm:ss.SS';
fprintf('\n elapsed time: %s\n', string(d));

fprintf('\nHawkes fit summary:\n');
fprintf('  NLL = %.2f,  BIC = %.2f\n', params.nll, params.BIC);
fprintf('  Baseline rates (events/s):\n');
for e=1:nE, fprintf('    %s: %.4f\n', eventnames{e}, params.mu(e)); end

if opt.Plot
    figure('Position',[100 100 1200 430]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(params.alpha); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('\alpha  (excitation amplitude)'); colormap(hot);

    nexttile;
    imagesc(1./params.beta); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('1/\beta  (decay timescale, s)'); colormap(hot);

    nexttile;
    br = params.branchingRatio;
    imagesc(br, [0 1]); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('\alpha/\beta  (branching ratio: fraction of triggered events)');
    colormap(gca, hot);

    sgtitle('Hawkes process fit  [row=triggered event, col=triggering event]');
end

end

function nll = hawkes_nll(x, eventTimes, nE, T)
% Negative log-likelihood for multivariate Hawkes process
    mu    = x(1:nE);
    alpha = reshape(x(nE+1:nE+nE^2),        nE, nE);
    beta  = reshape(x(nE+nE^2+1:nE+2*nE^2), nE, nE);

    loglik = 0;

    % fprintf('\n fitting hawkes nll ... ');
    % nbytes=  fprintf('%d/%d', 0, nE);



    for e = 1:nE
    % fprintf(repmat('\b',1,nbytes));
    % nbytes=  fprintf('%d/%d', e, nE);
        te = eventTimes{e};
        if isempty(te), continue; end
        % Integral term (compensator)
        compens = mu(e) * T;
        for f = 1:nE
            tf = eventTimes{f};
            if isempty(tf), continue; end
            % sum over triggering events: alpha/beta*(1-exp(-beta*(T-tf)))
            compens = compens + (alpha(e,f)/beta(e,f)) * ...
                      sum(1 - exp(-beta(e,f)*(T - tf(tf < T))));
        end
        % Sum log-intensity at each event time
        logIntens = 0;
        for k = 1:numel(te)
            lam = mu(e);
            for f = 1:nE
                tf_prev = eventTimes{f}(eventTimes{f} < te(k));
                if ~isempty(tf_prev)
                    lam = lam + alpha(e,f) * sum(exp(-beta(e,f)*(te(k)-tf_prev)));
                end
            end
            logIntens = logIntens + log(max(lam, 1e-10));
        end
        loglik = loglik + logIntens - compens;
    end
    nll = -loglik;
end


%% =========================================================================
function [simTimes] = pp_hawkes_simulate(params, eventnames, duration, varargin)
% PP_HAWKES_SIMULATE  Simulate from a fitted multivariate Hawkes process
%   Uses Ogata's thinning algorithm
%
%   [simTimes] = pp_hawkes_simulate(params, eventnames, duration)
%   params   : output of pp_hawkes_fit
%   duration : simulation duration in seconds
%
%   Optional name-value:
%     'Plot' : true/false (default true)

p = inputParser;
p.addParameter('Plot', true);
p.parse(varargin{:});
opt = p.Results;

nE    = numel(eventnames);
mu    = params.mu;
alpha = params.alpha;
beta  = params.beta;

simTimes = cell(nE, 1);
history  = cell(nE, 1);  % accumulated event times per type

t = 0;
while t < duration
    % Upper bound on total intensity (sum of all baselines + all past excitation)
    lambdaStar = sum(mu);
    for e = 1:nE
        for f = 1:nE
            if ~isempty(history{f})
                lambdaStar = lambdaStar + alpha(e,f) * ...
                    sum(exp(-beta(e,f)*(t - history{f})));
            end
        end
    end
    lambdaStar = max(lambdaStar, 1e-10);

    % Draw candidate next event time
    dt = -log(rand) / lambdaStar;
    t  = t + dt;
    if t > duration, break; end

    % Compute actual intensity at t
    lambdaActual = zeros(nE,1);
    for e = 1:nE
        lambdaActual(e) = mu(e);
        for f = 1:nE
            if ~isempty(history{f})
                lambdaActual(e) = lambdaActual(e) + alpha(e,f) * ...
                    sum(exp(-beta(e,f)*(t - history{f})));
            end
        end
    end
    lambdaTotal = sum(lambdaActual);

    % Accept/reject
    if rand <= lambdaTotal / lambdaStar
        % Which event type?
        eType = find(rand <= cumsum(lambdaActual/lambdaTotal), 1, 'first');
        history{eType}(end+1) = t;
        simTimes{eType}(end+1) = t;
    end
end

simTimes = cellfun(@(x) x(:), simTimes, 'UniformOutput', false);

if opt.Plot
    cols = lines(nE);
    figure('Position',[100 100 1400 500]);
    hold on;
    for e = 1:nE
        st = simTimes{e};
        plot(st, e*ones(size(st)), '|', 'Color', cols(e,:), ...
            'MarkerSize', 14, 'LineWidth', 1.5);
    end
    yticks(1:nE); yticklabels(eventnames);
    set(gca,'TickLabelInterpreter','none');
    xlabel('Time (s)');
    title(sprintf('Simulated Hawkes process  (duration=%.0fs)', duration));
    xlim([0 duration]); ylim([0.5 nE+0.5]);
    grid on;
end
end


%% =========================================================================
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


%% =========================================================================
function pp_summary(eventFrames, eventnames, stateMask, statenames, num_frames, fps)
% PP_SUMMARY  Run and display all point process analyses
%   Convenience wrapper — calls all functions with default parameters

fprintf('=== POINT PROCESS ANALYSIS SUITE ===\n\n');

fprintf('[1/5] Computing kernel-smoothed event rates...\n');
pp_rate(eventFrames, eventnames, num_frames, fps);

fprintf('[2/5] Computing inter-event interval distributions...\n');
pp_iei(eventFrames, eventnames, fps);

fprintf('[3/5] Computing cross-intensity functions...\n');
pp_crossintensity(eventFrames, eventnames, num_frames, fps);

fprintf('[4/5] Computing state-conditional event rates...\n');
pp_state_rates(eventFrames, eventnames, stateMask, statenames, fps);

fprintf('[5/5] Fitting point process GLM...\n');
pp_glm(eventFrames, eventnames, stateMask, statenames, num_frames, fps);

fprintf('\nDone. Note: pp_hawkes_fit and pp_hawkes_simulate must be run\n');
fprintf('separately as fitting can be slow for large datasets.\n');
end


%% =========================================================================
%  SHARED UTILITY: red-blue diverging colormap
function cmap = redblue_cmap(n)
if nargin < 1, n = 256; end
r = [linspace(0,1,n/2), ones(1,n/2)];
b = [ones(1,n/2), linspace(1,0,n/2)];
g = [linspace(0,1,n/2), linspace(1,0,n/2)];
cmap = [r' g' b'];
end
