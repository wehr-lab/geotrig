% pp_analysis_by_condition


statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
eventnames = {'failed_approach','contact_loss', ...
    'contact_gain','intercept','cricket_jump','rangemin'};

nStates = numel(statenames);
nEvents = numel(eventnames);

eventFrames={failed_approach_event_frames,contact_loss_event_frames, ...
    contact_gain_event_frames,intercept_event_frames,cricket_jump_event_frames,rangemin_event_frames};

fps=200;

%pp_summary(eventFrames, eventnames, stateMask, statenames, num_frames, fps)
clear condition_name

for c=1:4

    if c==1
        condition=laseron & dark; condition_name{c}='dark laser on';
    elseif c==2
        condition=laseron & light; condition_name{c}='light laser on';
    elseif c==3
        condition=~laseron & dark; condition_name{c}='dark laser off';
    elseif c==4
        condition=~laseron & light; condition_name{c}='light laser off';
    end


stateMask = [    condition & hotpursuit(:), condition & chase(:), condition & follow(:),...
     condition & stalk(:),  condition & wander(:),  condition & pause(:)];


    fprintf('Computing state-conditional event rates...\n');
    [stateRates, rateMat] = pp_state_rates(eventFrames, eventnames, stateMask, statenames, fps, 'title', condition_name{c});
    sgtitle(condition_name{c})
    rateMats{c}=rateMat;
end



figure('Position',[100 100 1200 860]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for c=1:4
    nexttile;
    imagesc(rateMats{c});
    colorbar; colormap(hot);
    set(gca,'XTick',1:nEvents,'XTickLabel',eventnames,'XTickLabelRotation',40,...
        'YTick',1:nStates,'YTickLabel',statenames,'TickLabelInterpreter','none');
    title(['Event rate (events/s) by state, ', condition_name{c}]);
end



figure('Position',[100 100 1200 860]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
A=[1 2 1 3];
B=[3 4 2 4];
for c=1:4
    a=A(c);
    b=B(c);
    nexttile;
    diffmat=rateMats{a}-rateMats{b};
    imagesc(diffmat);
    cmax=max(abs(diffmat(:)));
    clim(cmax*[-1 1])
    colorbar; colormap(hot);
    set(gca,'XTick',1:nEvents,'XTickLabel',eventnames,'XTickLabelRotation',40,...
        'YTick',1:nStates,'YTickLabel',statenames,'TickLabelInterpreter','none');
    title(['Event rate (events/s) by state, ', condition_name{a},' - ', condition_name{b}]);
end
colormap redblue_cmap(10)
