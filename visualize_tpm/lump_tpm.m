function [P_coarse, info] = lump_tpm(P_fine, partition, varargin)
% LUMP_TPM  Derive a coarse-grained transition probability matrix from a
% fine-grained TPM given a partition of states into clusters, and assess
% how well the lumping approximates a true Markov chain.
%
% USAGE:
%   [P_coarse, info] = lump_tpm(P_fine, partition)
%   [P_coarse, info] = lump_tpm(P_fine, partition, 'pi', pi_fine)
%
% INPUTS:
%   P_fine    - N x N transition probability matrix (rows sum to 1)
%   partition - N x 1 vector of cluster labels (integers 1..K), assigning
%               each fine state to a coarse state. Need not be contiguous
%               or sorted.
%
% OPTIONAL NAME-VALUE ARGS:
%   'pi'      - 1 x N (or N x 1) stationary distribution of P_fine. If not
%               supplied, it is computed automatically from P_fine.
%
% OUTPUTS:
%   P_coarse  - K x K lumped transition probability matrix, computed as a
%               stationary-distribution-weighted aggregation:
%                   P_coarse(a,b) = sum_{i in a, j in b} pi(i)*P_fine(i,j)
%                                   / sum_{i in a} pi(i)
%               i.e. the probability of moving from cluster a to cluster b
%               given that the chain is currently in cluster a, averaged
%               over the within-cluster stationary distribution.
%
%   info      - struct with diagnostic fields:
%       .pi_fine        - stationary distribution of P_fine (1 x N)
%       .pi_coarse      - stationary distribution of P_coarse (1 x K),
%                          equal to summing pi_fine within each cluster
%       .M              - N x K membership (indicator) matrix
%       .lumpability_err- K x K matrix giving, for each (a,b), the
%                          weighted standard deviation across i in cluster
%                          a of P(i -> cluster b). This is the classic
%                          "exact lumpability" check (Kemeny-Snell): if
%                          for every i in cluster a, sum_{j in b} P(i,j)
%                          is the same value, the lumping is exact and
%                          this entry is 0. Larger values indicate the
%                          lumped chain is a worse Markov approximation
%                          for that (a,b) pair.
%       .max_err        - scalar, max entry of lumpability_err (quick
%                          single-number summary of lumping quality)
%       .row_probs      - cell array, row_probs{a}{b} gives the N_a x 1
%                          vector of sum_{j in b} P(i,j) for each i in
%                          cluster a (the raw values info.lumpability_err
%                          is derived from -- useful for histograms)
%
% EXAMPLE:
%   P = [0.80 0.15 0.03 0.02;
%        0.10 0.85 0.03 0.02;
%        0.02 0.03 0.80 0.15;
%        0.02 0.03 0.10 0.85];
%   partition = [1 1 2 2];   % lump {1,2} -> coarse state 1, {3,4} -> coarse state 2
%   [Pc, info] = lump_tpm(P, partition);
%   disp(Pc)
%   disp(info.max_err)

    p = inputParser;
    addRequired(p, 'P_fine', @(x) ismatrix(x) && size(x,1)==size(x,2));
    addRequired(p, 'partition', @(x) isvector(x));
    addParameter(p, 'pi', [], @(x) isempty(x) || isvector(x));
    parse(p, P_fine, partition, varargin{:});

    P_fine = p.Results.P_fine;
    partition = p.Results.partition(:); % column vector
    N = size(P_fine,1);

    if numel(partition) ~= N
        error('lump_tpm:badPartition', ...
            'partition must have one entry per state (length %d), got %d.', N, numel(partition));
    end

    labels = unique(partition);
    K = numel(labels);

    % Remap arbitrary labels to 1..K for indexing convenience
    [~, ~, partIdx] = unique(partition);

    % --- Membership matrix M (N x K), one-hot ---
    M = zeros(N, K);
    for i = 1:N
        M(i, partIdx(i)) = 1;
    end

    % --- Stationary distribution of P_fine ---
    if isempty(p.Results.pi)
        pi_fine = stationary_distribution(P_fine);
    else
        pi_fine = p.Results.pi(:)'; % row vector
        if abs(sum(pi_fine) - 1) > 1e-6
            warning('lump_tpm:piNotNormalized', ...
                'Provided pi does not sum to 1 (sum = %.6f); normalizing.', sum(pi_fine));
            pi_fine = pi_fine / sum(pi_fine);
        end
    end

    % --- Coarse stationary distribution: sum pi within each cluster ---
    pi_coarse = pi_fine * M; % 1 x K

    % --- Coarse TPM via stationary-weighted aggregation ---
    % Numerator: pi_fine(i) * P_fine(i,j), aggregated to (cluster_a, cluster_b)
    weighted_P = (pi_fine(:) .* P_fine); % N x N, row i scaled by pi_fine(i)
    numer = M' * weighted_P * M;          % K x K
    denom = pi_coarse(:);                 % K x 1 (sum of pi over cluster a)

    P_coarse = numer ./ denom;            % broadcast K x K ./ K x 1

    % Guard against rows with zero stationary mass (shouldn't normally happen)
    zero_rows = denom < eps;
    if any(zero_rows)
        warning('lump_tpm:zeroMassCluster', ...
            'Cluster(s) %s have ~zero stationary probability; corresponding P_coarse rows set to uniform.', ...
            mat2str(find(zero_rows)'));
        P_coarse(zero_rows, :) = 1/K;
    end

    % Renormalize rows defensively (should already sum to ~1)
    P_coarse = P_coarse ./ sum(P_coarse, 2);

    % --- Lumpability error diagnostic ---
    % For each cluster a and cluster b, look at the per-state quantity
    %   r_i = sum_{j in b} P_fine(i,j),  for each i in cluster a
    % Exact lumpability requires r_i to be constant across i in a.
    % We report the stationary-weighted std of r_i as the error measure.
    lumpability_err = zeros(K,K);
    row_probs = cell(K,1);

    % Precompute, for every state i, its probability of landing in each cluster b
    P_to_cluster = P_fine * M; % N x K : P_to_cluster(i,b) = sum_{j in b} P_fine(i,j)

    for a = 1:K
        idx_a = find(partIdx == a);
        w = pi_fine(idx_a);
        w = w / sum(w); % weights within cluster a (stationary-conditional)

        row_probs{a} = cell(K,1);
        for b = 1:K
            r = P_to_cluster(idx_a, b); % N_a x 1
            row_probs{a}{b} = r;

            if numel(r) <= 1
                lumpability_err(a,b) = 0;
            else
                wmean = sum(w(:) .* r);
                wvar = sum(w(:) .* (r - wmean).^2);
                lumpability_err(a,b) = sqrt(wvar);
            end
        end
    end

    info = struct();
    info.pi_fine = pi_fine;
    info.pi_coarse = pi_coarse;
    info.M = M;
    info.lumpability_err = lumpability_err;
    info.max_err = max(lumpability_err(:));
    info.row_probs = row_probs;
    info.partIdx = partIdx; % normalized 1..K labels, in case original labels weren't

end


function pi_dist = stationary_distribution(P)
% STATIONARY_DISTRIBUTION  Compute the stationary distribution of an
% irreducible, aperiodic Markov chain with transition matrix P (rows sum
% to 1), as the left eigenvector for eigenvalue 1, normalized to sum to 1.

    N = size(P,1);
    [V, D] = eig(P');
    [~, idx] = min(abs(diag(D) - 1));
    pi_dist = real(V(:,idx))';
    pi_dist = pi_dist / sum(pi_dist);

    % Guard against sign/numerical issues
    if any(pi_dist < -1e-8)
        warning('stationary_distribution:negativeEntries', ...
            'Stationary distribution has small negative entries; clipping to 0 and renormalizing.');
        pi_dist(pi_dist < 0) = 0;
        pi_dist = pi_dist / sum(pi_dist);
    end
end
