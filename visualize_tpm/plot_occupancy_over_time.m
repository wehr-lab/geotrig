function plot_occupancy_over_time(P, T, varargin)
% PLOT_OCCUPANCY_OVER_TIME  Plot how probability mass evolves across
% states over time, by propagating an initial distribution forward under
% P (i.e. pi_t = pi_0 * P^t). Optionally group/color fine states by a
% partition to see how clusters behave individually and in aggregate.
%
% USAGE:
%   plot_occupancy_over_time(P, T)
%   plot_occupancy_over_time(P, T, 'pi0', pi0)
%   plot_occupancy_over_time(P, T, 'partition', partition)
%   plot_occupancy_over_time(P, T, 'partition', partition, 'pi0', pi0, 'stateNames', names)
%
% INPUTS:
%   P  - N x N transition probability matrix (rows sum to 1)
%   T  - number of time steps to propagate (t = 0, 1, ..., T)
%
% OPTIONAL NAME-VALUE ARGS:
%   'pi0'        - 1 x N initial distribution (default: point mass on
%                  state 1, i.e. pi0 = [1 0 0 ... 0])
%   'partition'  - N x 1 vector of cluster labels (1..K). If provided,
%                  produces TWO subplots:
%                     (top)    fine-state occupancy, lines colored by
%                              cluster membership
%                     (bottom) coarse-state occupancy (sum of fine states
%                              within each cluster), one line per cluster
%                  If omitted, produces a single plot of fine-state
%                  occupancy with default coloring.
%   'stateNames' - cell array of N strings for fine-state labels
%   'clusterNames' - cell array of K strings for coarse-cluster labels
%                  (only used if 'partition' is given)
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];
%   plot_occupancy_over_time(P, 30, 'partition', partition, ...
%       'stateNames', {'A1','A2','B1','B2'}, ...
%       'clusterNames', {'Cluster A','Cluster B'});

    p = inputParser;
    addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addRequired(p, 'T', @(x) isscalar(x) && x >= 1);
    addParameter(p, 'pi0', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'partition', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'stateNames', {}, @iscell);
    addParameter(p, 'clusterNames', {}, @iscell);
    parse(p, P, T, varargin{:});

    P = p.Results.P;
    T = p.Results.T;
    N = size(P,1);

    fps=200;

    if isempty(p.Results.pi0)
        pi0 = zeros(1,N);
        pi0(1) = 1;
    else
        pi0 = p.Results.pi0(:)';
        if abs(sum(pi0) - 1) > 1e-6
            warning('plot_occupancy_over_time:pi0NotNormalized', ...
                'pi0 does not sum to 1 (sum = %.6f); normalizing.', sum(pi0));
            pi0 = pi0 / sum(pi0);
        end
    end

    if isempty(p.Results.stateNames)
        stateNames = arrayfun(@(i) sprintf('S%d', i), 1:N, 'UniformOutput', false);
    else
        stateNames = p.Results.stateNames;
    end

    % --- Propagate distribution forward: pi_history(t+1,:) = pi0 * P^t ---
    pi_history = zeros(T+1, N);
    pi_history(1,:) = pi0;
    for t = 1:T
        pi_history(t+1,:) = pi_history(t,:) * P;
    end
    timeAxis = 0:T;
    timeAxis=timeAxis/fps;

    partition = p.Results.partition;
    cmap=flipud(parula( N+4));

    if isempty(partition)
        % --- Single plot: fine-state occupancy ---
        figure('Color','w');
        h=plot(timeAxis, pi_history, 'LineWidth', 2);
        for i=1:3
            set(h(i), 'color', cmap(i+2,:))
        end
        for i=4:6
            set(h(i), 'color', cmap(i+2,:))
        end
        xlabel('Time, s');
        ylabel('P(X_t = state)');
        title(sprintf('State Occupancy Over Time, starting from  %s', stateNames{find(pi0)}), 'Interpreter','none','FontSize',9);
        legend(stateNames, 'Location', 'eastoutside');
        grid on;
        ylim([0 1]);
        set(gcf, 'pos', [  1001         944         986         295])
        set(gca, 'fontsiz', 18)
    else
        % --- Two-panel plot: fine-state (colored by cluster) + coarse ---
        partition = partition(:);
        if numel(partition) ~= N
            error('plot_occupancy_over_time:badPartition', ...
                'partition must have one entry per state (length %d), got %d.', N, numel(partition));
        end
        [~, ~, partIdx] = unique(partition);
        K = max(partIdx);

        if isempty(p.Results.clusterNames)
            clusterNames = arrayfun(@(k) sprintf('Cluster %d', k), 1:K, 'UniformOutput', false);
        else
            clusterNames = p.Results.clusterNames;
        end

        % Coarse occupancy: sum fine-state probabilities within each cluster
        pi_coarse_history = zeros(T+1, K);
        for k = 1:K
            pi_coarse_history(:,k) = sum(pi_history(:, partIdx==k), 2);
        end

        % Color map: one base color per cluster, fine states within a
        % cluster get shades of that color
        clusterColors = lines(K);

        figure('Color','w');

        % --- Top: fine states, colored/shaded by cluster ---
        subplot(2,1,1);
        hold on;
        legendHandles = gobjects(N,1);
        for k = 1:K
            idx = find(partIdx == k);
            nInCluster = numel(idx);
            for ii = 1:nInCluster
                i = idx(ii);
                % Vary shade within cluster: lighter for later members
                shadeFrac = 0.3 + 0.5 * (ii-1) / max(nInCluster-1,1); % 0.3..0.8
                col = clusterColors(k,:) * shadeFrac + (1-shadeFrac)*[1 1 1];
                legendHandles(i) = plot(timeAxis, pi_history(:,i), '-', ...
                    'Color', col, 'LineWidth', 1.5);
            end
        end
        xlabel('Time step t');
        ylabel('P(X_t = state)');
        title('Fine-State Occupancy (colored by cluster)');
        legend(legendHandles, stateNames, 'Location', 'eastoutside');
        grid on;
        ylim([0 1]);
        hold off;

        % --- Bottom: coarse (lumped) states ---
        subplot(2,1,2);
        hold on;
        for k = 1:K
            plot(timeAxis, pi_coarse_history(:,k), '-', ...
                'Color', clusterColors(k,:), 'LineWidth', 2.5);
        end
        xlabel('Time step t');
        ylabel('P(X_t \in cluster)');
        title('Coarse (Lumped) Cluster Occupancy');
        legend(clusterNames, 'Location', 'eastoutside');
        grid on;
        ylim([0 1]);
        hold off;
    end

end

