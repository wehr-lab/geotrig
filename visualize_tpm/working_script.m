% visualize TPM working script

% P = [0.80 0.15 0.03 0.02;
%      0.10 0.85 0.03 0.02;
%      0.02 0.03 0.80 0.15;
%      0.02 0.03 0.10 0.85];

P=TPMdark;
%TPM was computed by PlotTPM, using 4 chase states and all the rest unsplit
statenames'

% partition = [1 2 2 3 3 4 5 6 7 8 9 10 ];
%partition =  [1 2 2 2 2 3 4 5 6 7 8 9  ];
% partition =   [1 1 1 1 1 1 1 3 2 2 3 2 ];
partition = [ 1 1 1 2 2 2];
partition = [ 1 1 1 1 2 2];

pi0=contains(statenames, 'pause'); %0s for all states except 1 at the match

clusterNames = derive_cluster_names(statenames, partition)

[Pc, info] = lump_tpm(P, partition);
disp(Pc); disp(info.max_err)

traj_fine = simulate_markov_chain(P, 1000);
traj_coarse = relabel_trajectory(traj_fine, partition);

plot_tpm_circle_grouped(P, statenames, partition, ...
    struct('clusterNames', {clusterNames}));

plot_tpm_circle_grouped(P, statenames);

plot_occupancy_over_time(P, 100, 'partition', partition, ...
    'stateNames', statenames, 'clusterNames', clusterNames);

plot_occupancy_over_time(P, 200,    'stateNames', statenames, 'pi0', pi0);

[Pc, info] = lump_tpm(P, partition);
fprintf('Max lumpability error: %.4f\n', info.max_err);

plot_mfpt_heatmap(P, statenames, 'partition', partition);

for c=1:4
    plot_mfpt_heatmap(results(c).TPM, statenames, 'title', results(c).condition_name);
   % plot_mfpt_heatmap(results(c).TPM, statenames, 'title', results(c).condition_name, 'partition', partition);
end

%% test clusters using MFPT
c=3;
MFPT=plot_mfpt_heatmap(results(c).TPM, statenames, 'title', results(c).condition_name);

pursuitidx=find(partition==1)
searchidx=find(partition==2)
find(partition==1)
pursuit=MFPT(pursuitidx, pursuitidx);
fprintf('\nmean(pursuit(:)) %.f', mean(pursuit(:)))
search=MFPT(searchidx, searchidx);
fprintf('\nmean(search(:)) %.f', mean(search(:)))
inter1=MFPT(pursuitidx, searchidx);
inter2=MFPT(searchidx, pursuitidx);
fprintf('\ninter-1 %.f', mean(inter1(:)))
fprintf('\ninter-2 %.f', mean(inter2(:)))
fprintf('\nmean inter/intra ratio %.2f', mean([inter1(:); inter2(:)])/mean([search(:); pursuit(:)]) )

%% 1. Spectrum suggests how many clusters to try
% plots eigenvalues of P both in the complex plane (with the unit circle)
% and as a sorted bar chart of magnitudes. The thing to look for is a gap
% in |λ_k| after the first few eigenvalues — if |λ_2|,...,|λ_K| are all
% close to 1 and then there's a drop to |λ_{K+1}|, that suggests K natural
% slow/metastable clusters, giving you a principled starting guess for how
% many coarse states to aim for.
c=3;
eigvals=plot_tpm_spectrum(results(c).TPM);



%pi0=zeros(size(statenames));

plot_tpm_sankey(P, 6, 'partition', partition, 'stateNames', statenames, 'pi0', pi0);

plot_tpm_sankey(P, 6,  'stateNames', statenames, 'pi0', pi0);

traj_all = plot_sample_trajectories(P, 200, 30, 'partition', partition, ...
    'stateNames', statenames, 'clusterNames', clusterNames, 'x0', 1);

traj_all = plot_sample_trajectories(P, 2000, 5,  'stateNames', statenames);



traj = simulate_markov_chain(P, 50000);
i=1;
rt_state1  = plot_return_time_distribution(P, i, 'traj', traj, 'label', statenames{i});
rt_clusterA = plot_return_time_distribution(P, i, 'partition', partition, 'traj', traj, 'label', clusterNames{i});














