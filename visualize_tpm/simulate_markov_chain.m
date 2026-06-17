function traj = simulate_markov_chain(P, T, x0)
% SIMULATE_MARKOV_CHAIN  Simulate a discrete-time Markov chain trajectory.
%
% USAGE:
%   traj = simulate_markov_chain(P, T)
%   traj = simulate_markov_chain(P, T, x0)
%
% INPUTS:
%   P  - N x N transition probability matrix (rows sum to 1)
%   T  - number of time steps to simulate (length of trajectory)
%   x0 - (optional) initial state (1..N). If omitted, drawn from the
%        stationary distribution of P.
%
% OUTPUT:
%   traj - 1 x T vector of state indices (1..N)
%
% NOTE:
%   Because this works directly from P_fine, the resulting trajectory can
%   be reused for any candidate lumping via relabel_trajectory, without
%   re-simulating.

    N = size(P,1);

    if nargin < 3 || isempty(x0)
        % Lazy stationary distribution via power iteration (avoids needing
        % lump_tpm.m's eig-based version as a dependency)
        pi0 = ones(1,N)/N;
        for k = 1:1000
            pi1 = pi0 * P;
            if norm(pi1 - pi0, 1) < 1e-12
                pi0 = pi1;
                break;
            end
            pi0 = pi1;
        end
        x0 = randsample_local(1:N, pi0);
    end

    traj = zeros(1,T);
    traj(1) = x0;
    cumP = cumsum(P, 2); % precompute cumulative rows for fast sampling

    for t = 2:T
        prev = traj(t-1);
        r = rand();
        traj(t) = find(cumP(prev,:) >= r, 1, 'first');
    end
end


function s = randsample_local(vals, weights)
% Minimal weighted sampling without requiring Statistics Toolbox
    weights = weights / sum(weights);
    c = cumsum(weights);
    r = rand();
    idx = find(c >= r, 1, 'first');
    s = vals(idx);
end
