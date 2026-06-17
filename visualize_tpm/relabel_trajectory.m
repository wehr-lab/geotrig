function traj_coarse = relabel_trajectory(traj_fine, partition)
% RELABEL_TRAJECTORY  Map a fine-grained state sequence to a coarse-grained
% sequence given a partition, so simulated trajectories can be reused
% across different candidate lumpings without re-simulating.
%
% USAGE:
%   traj_coarse = relabel_trajectory(traj_fine, partition)
%
% INPUTS:
%   traj_fine  - 1 x T (or T x 1) vector of fine-state indices (1..N),
%                e.g. output of simulate_markov_chain
%   partition  - N x 1 vector of cluster labels (integers 1..K), same
%                convention as used in lump_tpm. Arbitrary label values
%                are remapped to 1..K.
%
% OUTPUT:
%   traj_coarse - same shape as traj_fine, with each fine-state index
%                  replaced by its coarse cluster index (1..K)
%
% EXAMPLE:
%   traj_fine = [1 2 1 3 4 3 4 4 1];
%   partition = [1 1 2 2]; % states {1,2}->cluster 1, {3,4}->cluster 2
%   traj_coarse = relabel_trajectory(traj_fine, partition);
%   % traj_coarse = [1 1 1 2 2 2 2 2 1]

    partition = partition(:);
    [~, ~, partIdx] = unique(partition); % normalize labels to 1..K

    if any(traj_fine(:) < 1) || any(traj_fine(:) > numel(partition))
        error('relabel_trajectory:badIndex', ...
            'traj_fine contains indices outside the range 1..%d implied by partition.', numel(partition));
    end

    traj_coarse = reshape(partIdx(traj_fine(:)), size(traj_fine));
end
