function [H_mean, H_std, H_vals] = hopkins_statistic(X, m, n_iter, metric)
% HOPKINS_STATISTIC  Compute the Hopkins statistic for a dataset X.
%
% Usage:
%   [H_mean, H_std, H_vals] = hopkins_statistic(X)
%   [H_mean, H_std, H_vals] = hopkins_statistic(X, m, n_iter, metric)
%
% Inputs:
%   X       - (n x d) matrix, one row per sample (e.g. rows of a TPM)
%   m       - number of points sampled per trial (default: max(2, floor(n/10)))
%   n_iter  - number of Monte-Carlo repetitions (default: 200)
%   metric  - 'euclidean' (default) or 'kl_sym'
%
% Outputs:
%   H_mean  - mean Hopkins statistic over n_iter trials
%   H_std   - standard deviation across trials
%   H_vals  - (n_iter x 1) vector of per-trial H values
%
% Interpretation:
%   H ~ 0.5        no cluster tendency (random spatial distribution)
%   0.5 < H < 0.75 weak tendency
%   0.75 <= H < 0.9 moderate cluster tendency
%   H >= 0.9       strong cluster tendency
%   H ~ 1.0        perfectly clustered
%
% Notes on applying to a TPM:
%   Pass the rows of your transition matrix as X.  Each row is a
%   probability distribution over states, so 'kl_sym' is the more
%   principled metric.  With only n=6 states variance is high; use
%   n_iter >= 200 and treat the result as indicative.

    [n, d] = size(X);

    if nargin < 2 || isempty(m),      m      = max(2, floor(n/10)); end
    if nargin < 3 || isempty(n_iter), n_iter = 200;                 end
    if nargin < 4 || isempty(metric), metric = 'euclidean';         end

    if m >= n
        error('hopkins_statistic: m (%d) must be < n (%d).', m, n);
    end

    % Bounding box for uniform random sampling
    lo   = min(X, [], 1);          % (1 x d)
    hi   = max(X, [], 1);
    span = hi - lo;
    span(span == 0) = 1;           % avoid zero-width dimensions

    dist_fn = pick_metric(metric);

    H_vals = nan(n_iter, 1);

    for iter = 1:n_iter
        % --- sample m real points without replacement ---
        idx     = randperm(n, m);
        X_samp  = X(idx, :);
        X_rest  = X(setdiff(1:n, idx), :);   % remaining n-m points

        % nearest-neighbour distances: real sample -> rest of X
        u = nn_dist(X_samp, X_rest, dist_fn);    % (m x 1)

        % --- sample m synthetic points uniformly in bounding box ---
        X_rand = lo + rand(m, d) .* span;

        % nearest-neighbour distances: synthetic -> all of X
        w = nn_dist(X_rand, X, dist_fn);          % (m x 1)

        denom = sum(u) + sum(w);
        if denom > 0
            H_vals(iter) = sum(w) / denom;
        end
    end

    H_vals = H_vals(~isnan(H_vals));
    H_mean = mean(H_vals);
    H_std  = std(H_vals);
end


% ── helpers ──────────────────────────────────────────────────────────────────

function dmin = nn_dist(A, B, dist_fn)
% For each row in A, find the distance to its nearest neighbour in B.
    dmin = zeros(size(A, 1), 1);
    for i = 1:size(A, 1)
        dists = dist_fn(A(i, :), B);   % (size(B,1) x 1)
        dmin(i) = min(dists);
    end
end


function fn = pick_metric(metric)
    switch lower(metric)
        case 'euclidean'
            fn = @(row, mat) sqrt(sum((mat - row).^2, 2));

        case 'kl_sym'
            % Symmetric KL divergence — appropriate for rows of a TPM.
            % Adds a small epsilon to avoid log(0).
            eps = 1e-10;
            fn = @(row, mat) kl_sym_fn(row, mat, eps);

        otherwise
            error('hopkins_statistic: unknown metric "%s". Use "euclidean" or "kl_sym".', metric);
    end
end


function d = kl_sym_fn(row, mat, eps)
% Symmetric KL: D_sym(p||q) = D_KL(p||q) + D_KL(q||p)
    p = row + eps;                          % (1 x d)
    Q = mat + eps;                          % (m x d)
    kl_pq = sum(p .* (log(p) - log(Q)), 2);
    kl_qp = sum(Q .* (log(Q) - log(p)), 2);
    d = kl_pq + kl_qp;                     % (m x 1)
end