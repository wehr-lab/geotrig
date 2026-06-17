function returnTimes = plot_return_time_distribution(P, targetStates, varargin)
% PLOT_RETURN_TIME_DISTRIBUTION  Estimate and plot the distribution of
% return times to one or more states (or clusters of states) via
% simulation: starting a long trajectory, every time the chain visits
% (one of) the target state(s), record the number of steps since the
% previous visit.
%
% USAGE:
%   returnTimes = plot_return_time_distribution(P, targetStates)
%   returnTimes = plot_return_time_distribution(P, targetStates, 'T', 20000)
%   returnTimes = plot_return_time_distribution(P, partition, 'partition', partition, ...)
%
% INPUTS:
%   P            - N x N transition probability matrix (rows sum to 1)
%   targetStates - EITHER:
%                     (a) a vector of one or more state indices (1..N)
%                         whose return-time distribution to compute, OR
%                     (b) a single CLUSTER LABEL (scalar) if 'partition'
%                         is also provided -- in which case "return" means
%                         "return to ANY state in this cluster" (using the
%                         lumped/macro view of recurrence).
%
% OPTIONAL NAME-VALUE ARGS:
%   'T'          - length of simulated trajectory (default 20000). Longer
%                  trajectories give smoother estimates, especially for
%                  rare/long return times.
%   'partition'  - N x 1 vector of cluster labels (1..K), using the same
%                  [~,~,partIdx]=unique(partition) convention as the
%                  other functions. If provided AND targetStates is a
%                  scalar, targetStates is interpreted as a cluster label
%                  (after the unique() remapping) and "return" is to any
%                  state in that cluster.
%   'traj'       - a pre-simulated trajectory (1 x T) to reuse instead of
%                  simulating a new one (e.g. from simulate_markov_chain),
%                  for consistency across multiple calls/partitions.
%   'maxLag'     - cap the x-axis / histogram at this many steps (default:
%                  99th percentile of observed return times, for display
%                  only -- all return times are still used for the mean).
%   'theoretical'- if true (default), overlay the theoretical mean return
%                  time 1/pi(target) as a vertical line, where pi is the
%                  stationary distribution and pi(target) is summed over
%                  the target state(s)/cluster. Requires 'partition' not
%                  required for this -- works for both single states and
%                  clusters.
%   'label'      - string used in the plot title/legend to describe the
%                  target (default: auto-generated from targetStates)
%
% OUTPUT:
%   returnTimes  - vector of observed return times (gaps between
%                  consecutive visits to the target set), in steps
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%
%   % Return time to state 1 alone:
%   rt1 = plot_return_time_distribution(P, 1, 'label', 'State A1');
%
%   % Return time to cluster {1,2} (using a partition):
%   partition = [1 1 2 2];
%   rtA = plot_return_time_distribution(P, 1, 'partition', partition, ...
%       'label', 'Cluster A');
%
%   % Reuse the same trajectory for both clusters:
%   traj = simulate_markov_chain(P, 50000);
%   rtA = plot_return_time_distribution(P, 1, 'partition', partition, 'traj', traj, 'label','Cluster A');
%   rtB = plot_return_time_distribution(P, 2, 'partition', partition, 'traj', traj, 'label','Cluster B');

    p = inputParser;
    addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addRequired(p, 'targetStates', @(x) isvector(x));
    addParameter(p, 'T', 20000, @(x) isscalar(x) && x>=10);
    addParameter(p, 'partition', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'traj', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'maxLag', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'theoretical', true, @islogical);
    addParameter(p, 'label', '', @ischar);
    parse(p, P, targetStates, varargin{:});

    P = p.Results.P;
    targetStates = p.Results.targetStates;
    N = size(P,1);
    partition = p.Results.partition;

    % --- Resolve target set (which fine states count as "the target") ---
    if ~isempty(partition)
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_return_time_distribution:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);

        if isscalar(targetStates)
            clusterLabel = targetStates;
            targetSet = find(partIdx == clusterLabel);
            if isempty(targetSet)
                error('plot_return_time_distribution:badCluster', ...
                    'Cluster label %d not found (valid range 1..%d after unique()).', clusterLabel, max(partIdx));
            end
        else
            % Vector of fine states given alongside a partition: treat as
            % an explicit (possibly cross-cluster) target set
            targetSet = targetStates(:)';
        end
    else
        targetSet = targetStates(:)';
    end

    if any(targetSet < 1) || any(targetSet > N)
        error('plot_return_time_distribution:badTarget', ...
            'Target state indices must be in 1..%d.', N);
    end

    % --- Get trajectory (simulate or reuse) ---
    if ~isempty(p.Results.traj)
        traj = p.Results.traj;
    else
        traj = simulate_markov_chain(P, p.Results.T);
    end

    % --- Compute return times: gaps between consecutive visits to targetSet ---
    inTarget = ismember(traj, targetSet);
    visitTimes = find(inTarget);

    if numel(visitTimes) < 2
        warning('plot_return_time_distribution:tooFewVisits', ...
            'Target set was visited %d time(s) in the trajectory; cannot estimate a distribution. Try increasing T.', numel(visitTimes));
        returnTimes = [];
        return;
    end

    returnTimes = diff(visitTimes); % gaps in steps between consecutive visits

    % --- Label ---
    label = p.Results.label;
    if isempty(label)
        if isscalar(targetStates) && ~isempty(partition)
            label = sprintf('Cluster %d', targetStates);
        else
            label = sprintf('State(s) %s', mat2str(targetSet));
        end
    end

    % --- Theoretical mean return time via stationary distribution ---
    if p.Results.theoretical
        pi_fine = stationary_distribution_local(P);
        pi_target = sum(pi_fine(targetSet));
        meanReturnTheory = 1 / pi_target;
    end

    % --- Plot ---
    figure('Color','w');

    maxLag = p.Results.maxLag;
    if isempty(maxLag)
        maxLag = prctile(returnTimes, 99);
    end

    edges = 0.5:1:(maxLag+0.5);
    histogram(returnTimes, edges, 'Normalization','pdf', ...
        'FaceColor', [0.2 0.4 0.8], 'EdgeColor','none');
    hold on;

    meanEmpirical = mean(returnTimes);
    yl = ylim;
    plot([meanEmpirical meanEmpirical], yl, 'b-', 'LineWidth', 2);
    txt1 = sprintf('Empirical mean = %.2f', meanEmpirical);

    if p.Results.theoretical
        plot([meanReturnTheory meanReturnTheory], yl, 'r--', 'LineWidth', 2);
        txt2 = sprintf('Theoretical mean (1/\\pi) = %.2f', meanReturnTheory);
        legend({'Return time distribution', txt1, txt2}, 'Location','best');
    else
        legend({'Return time distribution', txt1}, 'Location','best');
    end

    xlabel('Return time (steps)');
    ylabel('Probability density');
    title(sprintf('Return Time Distribution: %s (n = %d visits)', label, numel(visitTimes)));
    xlim([0 maxLag]);
    grid on;
    hold off;
end


function pi_dist = stationary_distribution_local(P)
% Local copy of stationary distribution computation (eig-based), so this
% function has no dependency on lump_tpm.m.
    N = size(P,1);
    [V, D] = eig(P');
    [~, idx] = min(abs(diag(D) - 1));
    pi_dist = real(V(:,idx))';
    pi_dist = pi_dist / sum(pi_dist);

    if any(pi_dist < -1e-8)
        pi_dist(pi_dist < 0) = 0;
        pi_dist = pi_dist / sum(pi_dist);
    end
end
