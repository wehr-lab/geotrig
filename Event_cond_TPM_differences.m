%  Event_cond_TPM_differences
% script to compare event-conditioned TPMs across conditions


% if     c==1; condition = laseron &  dark;  condition_name = 'dark laser on';
% elseif c==2; condition = laseron & ~dark;  condition_name = 'light laser on';
% elseif c==3; condition = ~laseron &  dark; condition_name = 'dark laser off';
% elseif c==4; condition = ~laseron & ~dark; condition_name = 'light laser off';
% else

c=1;
Event_conditioned_TPMs_NOTA
results1=results;
c=2;
Event_conditioned_TPMs_NOTA
results2=results;
c=3;
Event_conditioned_TPMs_NOTA
results3=results;
c=4;
Event_conditioned_TPMs_NOTA
results4=results;

%% Circle diagrams

%dark - light (laser off)
figure('Position',[100 100 1600 1000]);
for e = 1:nEvents
    subplot(2,3,e)

    diffMatdarkon = results1(e).TPM_post - results1(e).TPM_baseline_excl;
    diffMatdarkoff = results3(e).TPM_post - results3(e).TPM_baseline_excl;
    diffMatlighton = results2(e).TPM_post - results2(e).TPM_baseline_excl;
    diffMatlightoff = results4(e).TPM_post - results4(e).TPM_baseline_excl;
    diffMat=diffMatdarkoff-diffMatlightoff;

    plot_tpm_diff_circle(diffMat, statenames, ...
        'Title', sprintf('%s (n=%d)', eventnames{e}, results(e).n_post_transitions), ...
        'MinAbsDiff', 0.001);
end
sgtitle(sprintf('TPM_{post} - TPM_{base}  [%s]', ...
    'dark-light (laser off)'));

%laseron - laseroff (dark)
figure('Position',[100 100 1600 1000]);
for e = 1:nEvents
    subplot(2,3,e)

    diffMatdarkon = results1(e).TPM_post - results1(e).TPM_baseline_excl;
    diffMatdarkoff = results3(e).TPM_post - results3(e).TPM_baseline_excl;
    diffMatlighton = results2(e).TPM_post - results2(e).TPM_baseline_excl;
    diffMatlightoff = results4(e).TPM_post - results4(e).TPM_baseline_excl;
    diffMat=diffMatdarkon - diffMatdarkoff;

    plot_tpm_diff_circle(diffMat, statenames, ...
        'Title', sprintf('%s (n=%d)', eventnames{e}, results(e).n_post_transitions), ...
        'MinAbsDiff', 0.001);
end
sgtitle(sprintf('TPM_{post} - TPM_{base}  [%s]', ...
    'laseron - laseroff (dark)'));

