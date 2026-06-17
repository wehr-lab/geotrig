function clusterNames = derive_cluster_names(stateNames, partition)
% DERIVE_CLUSTER_NAMES  Build cluster names from fine-state names and a
% partition vector, avoiding manual mismatches between cluster order and
% partition labels.
%
% USAGE:
%   clusterNames = derive_cluster_names(stateNames, partition)
%
% INPUTS:
%   stateNames - cell array of N strings (fine-state names)
%   partition  - N x 1 vector of cluster labels (arbitrary integers,
%                not necessarily 1..K or contiguous)
%
% OUTPUT:
%   clusterNames - K x 1 cell array, where clusterNames{k} corresponds to
%                  the cluster with label labels(k), labels = unique(partition).
%                  Each name is the member state names joined with '+',
%                  e.g. {'A1','A2'} -> 'A1+A2'.
%
% NOTE:
%   The ordering of clusterNames matches MATLAB's `unique`, i.e. cluster k
%   corresponds to the k-th smallest distinct value in `partition`. This
%   is the SAME ordering that lump_tpm, relabel_trajectory, and the
%   plotting functions use internally (they all call
%   [~, ~, partIdx] = unique(partition) ), so clusterNames{k} will line up
%   correctly with coarse-state index k everywhere.
%
% EXAMPLE:
%   stateNames = {'A1','A2','B1','B2'};
%   partition  = [1 1 2 2];
%   clusterNames = derive_cluster_names(stateNames, partition);
%   % clusterNames = {'A1+A2'; 'B1+B2'}
%
%   % Works even if labels are out of order / non-contiguous:
%   partition2 = [5 5 1 1];
%   clusterNames2 = derive_cluster_names(stateNames, partition2);
%   % labels = unique(partition2) = [1 5], so:
%   % clusterNames2 = {'B1+B2'; 'A1+A2'}   (cluster 1 = label 1 = {B1,B2})

    partition = partition(:);
    if numel(partition) ~= numel(stateNames)
        error('derive_cluster_names:sizeMismatch', ...
            'stateNames (length %d) and partition (length %d) must have the same length.', ...
            numel(stateNames), numel(partition));
    end

    [~, ~, partIdx] = unique(partition);
    K = max(partIdx);

    clusterNames = cell(K,1);
    for k = 1:K
        members = stateNames(partIdx == k);
        clusterNames{k} = strjoin(members, '+');
    end
end
