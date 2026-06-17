function traj_all = plot_sample_trajectories(P, T, nTraj, varargin)
% PLOT_SAMPLE_TRAJECTORIES  Simulate and plot one or more sample
% trajectories from a Markov chain, as step plots of state index over
% time. Optionally overlay a coarse-grained view based on a partition.
%
% USAGE:
%   traj_all = plot_sample_trajectories(P, T, nTraj)
%   traj_all = plot_sample_trajectories(P, T, nTraj, 'partition', partition)
%   traj_all = plot_sample_trajectories(P, T, nTraj, 'x0', x0)
%
% INPUTS:
%   P     - N x N transition probability matrix (rows sum to 1)
%   T     - number of time steps to simulate
%   nTraj - number of independent trajectories to simulate and overlay
%
% OPTIONAL NAME-VALUE ARGS:
%   'x0'           - initial state (1..N) for ALL trajectories. If
%                    omitted, each trajectory starts from a state drawn
%                    from the stationary distribution (independently).
%   'partition'    - N x 1 vector of cluster labels (1..K), using the
%                    same [~,~,partIdx]=unique(partition) convention as
%                    the other functions. If provided, produces TWO
%                    stacked subplots:
%                       (top)    fine-state trajectories, y-axis ordered
%                                so states are grouped by cluster (with
%                                horizontal band shading per cluster)
%                       (bottom) coarse-state trajectories (each fine
%                                trajectory relabeled via
%                                relabel_trajectory)
%                    If omitted, a single plot of fine-state trajectories
%                    is produced.
%   'stateNames'   - cell array of N strings for fine-state labels
%   'clusterNames' - cell array of K strings for coarse-cluster labels
%                    (only used if 'partition' is given)
%   'colors'       - nTraj x 3 matrix of RGB colors, one per trajectory
%                    (default: lines(nTraj))
%   'reuseTraj'    - N_traj-length cell array of pre-simulated fine
%                    trajectories to plot instead of simulating new ones
%                    (e.g. from simulate_markov_chain). Useful for
%                    comparing partitions on the EXACT SAME trajectories.
%                    If provided, nTraj and x0 are ignored.
%
% OUTPUT:
%   traj_all - 1 x nTraj cell array of fine-state trajectories (each
%              1 x T), so you can reuse them later via 'reuseTraj' for a
%              different partition.
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];
%   names = {'A1','A2','B1','B2'};
%   clusterNames = derive_cluster_names(names, partition);
%
%   traj_all = plot_sample_trajectories(P, 200, 3, 'partition', partition, ...
%       'stateNames', names, 'clusterNames', clusterNames, 'x0', 1);
%
%   % Reuse the same trajectories for a different partition:
%   plot_sample_trajectories(P, 200, 3, 'reuseTraj', traj_all, ...
%       'partition', [1 2 2 2], 'stateNames', names, ...
%       'clusterNames', derive_cluster_names(names, [1 2 2 2]));

    p = inputParser;
    addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addRequired(p, 'T', @(x) isscalar(x) && x>=1);
    addRequired(p, 'nTraj', @(x) isscalar(x) && x>=1);
    addParameter(p, 'x0', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'partition', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'stateNames', {}, @iscell);
    addParameter(p, 'clusterNames', {}, @iscell);
    addParameter(p, 'colors', [], @(x) isempty(x) || (ismatrix(x) && size(x,2)==3));
    addParameter(p, 'reuseTraj', {}, @iscell);
    parse(p, P, T, nTraj, varargin{:});

    P = p.Results.P;
    T = p.Results.T;
    N = size(P,1);

    if isempty(p.Results.stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    else
        stateNames = p.Results.stateNames;
    end

    % --- Get trajectories (simulate or reuse) ---
    if ~isempty(p.Results.reuseTraj)
        traj_all = p.Results.reuseTraj;
        nTraj = numel(traj_all);
        T = numel(traj_all{1});
    else
        traj_all = cell(1, nTraj);
        for k = 1:nTraj
            if isempty(p.Results.x0)
                traj_all{k} = simulate_markov_chain(P, T);
            else
                traj_all{k} = simulate_markov_chain(P, T, p.Results.x0);
            end
        end
    end

    if isempty(p.Results.colors)
        colors = lines(nTraj);
    else
        colors = p.Results.colors;
    end

    timeAxis = 0:(T-1);
    partition = p.Results.partition;

    if isempty(partition)
        % --- Single plot: fine-state trajectories ---
        figure('Color','w');
        hold on;
        for k = 1:nTraj
            stairs(timeAxis, traj_all{k}, '-', 'Color', colors(k,:), 'LineWidth', 1.2);
            stairs(timeAxis(end)+(T/100)*randn, traj_all{k}(end), 'o', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:));
            stairs(timeAxis(1), traj_all{k}(1), 'o', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:));
            
        end
        yticks(1:N);
        yticklabels(stateNames);
        ylim([0.5, N+0.5]);
        xlabel('Time step t');
        ylabel('State');
        title('Sample Trajectories (Fine States)');
        grid on;
        if nTraj > 1
            legend(arrayfun(@(k) sprintf('Trajectory %d',k), 1:nTraj, 'UniformOutput', false), ...
                'Location','eastoutside');
        end
        hold off;
    else
        % --- Two-panel plot: fine (grouped/shaded) + coarse ---
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_sample_trajectories:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);
        K = max(partIdx);

        if isempty(p.Results.clusterNames)
            clusterNames = arrayfun(@(k) sprintf('Cluster %d', k), 1:K, 'UniformOutput', false);
        else
            clusterNames = p.Results.clusterNames;
        end

        % --- Reorder fine states so same-cluster states are adjacent on
        % the y-axis (cluster-major order), for cleaner band shading ---
        [~, order] = sort(partIdx); % order(ii) = original state index at display position ii
        displayPos = zeros(N,1);    % displayPos(i) = y-position of original state i
        displayPos(order) = 1:N;
        clusterOfPos = partIdx(order); % cluster label at each display position

        clusterColors = lines(K);

        figure('Color','w');

        % --- Top: fine-state trajectories, y-axis grouped by cluster ---
        ax1 = subplot(2,1,1);
        hold on;

        % Shaded horizontal bands per cluster
        for pos = 1:N
            cl = clusterOfPos(pos);
            yBand = [pos-0.5, pos+0.5];
            xBand = [timeAxis(1), timeAxis(end)];
            patch([xBand(1) xBand(2) xBand(2) xBand(1)], ...
                  [yBand(1) yBand(1) yBand(2) yBand(2)], ...
                  clusterColors(cl,:), 'FaceAlpha', 0.08, 'EdgeColor','none');
        end

        for k = 1:nTraj
            traj_pos = displayPos(traj_all{k}); % map fine-state index -> display position
            stairs(timeAxis, traj_pos, '-', 'Color', colors(k,:), 'LineWidth', 1.2);
        end

        yticks(1:N);
        yticklabels(stateNames(order));
        ylim([0.5, N+0.5]);
        xlabel('Time step t');
        ylabel('State (grouped by cluster)');
        title('Sample Trajectories (Fine States, grouped by cluster)');
        grid on;
        hold off;

        % --- Bottom: coarse-state trajectories ---
        ax2 = subplot(2,1,2);
        hold on;

        for k = 1:nTraj
            traj_coarse = relabel_trajectory(traj_all{k}, partition);
            stairs(timeAxis, traj_coarse, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
        end

        yticks(1:K);
        yticklabels(clusterNames);
        ylim([0.5, K+0.5]);
        xlabel('Time step t');
        ylabel('Cluster');
        title('Sample Trajectories (Coarse Clusters)');
        grid on;
        hold off;

        if nTraj > 1
            legend(ax2, arrayfun(@(k) sprintf('Trajectory %d',k), 1:nTraj, 'UniformOutput', false), ...
                'Location','eastoutside');
        end

        linkaxes([ax1 ax2], 'x');
    end
end
