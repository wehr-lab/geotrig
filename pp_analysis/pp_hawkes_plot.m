function  pp_hawkes_plot( eventnames, params, varargin)
%just the plot from hawkes_fit without the fitting
if nargin==3 condition_name=varargin{1}; else condition_name=''; end
nE = numel(eventnames);
    figure('Position',[100 100 1300 430]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(params.alpha); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title(sprintf('\alpha  (excitation amplitude) %s', condition_name));
    colormap(hot);

    nexttile;
    imagesc(1./params.beta); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title('1/\beta  (decay timescale, s)'); 
    colormap(hot);

    nexttile;
    br = params.branchingRatio;
    imagesc(br, [0 1]); colorbar;
    set(gca,'XTick',1:nE,'XTickLabel',eventnames,'XTickLabelRotation',40,...
            'YTick',1:nE,'YTickLabel',eventnames,'TickLabelInterpreter','none');
    title(sprintf('\alpha/\beta  (branching ratio: fraction of triggered events)\n%s', condition_name));
    colormap(gca, hot);
    xlabel('triggering event')
    ylabel('received input')
    

    rho=max(abs(eig(params.branchingRatio)));
    sgtitle(sprintf('Hawkes process fit  [row=triggered event, col=triggering event], rho =%.2f', rho));
fprintf('\ncriticality: rho =%.4f', rho)
fprintf(['\nρ(K) < 1 → the process is stationary; cascades die out', ...
'\nρ(K) = 1 → critical; cascades can persist indefinitely (unit root behavior)', ...
'\nρ(K) > 1 → explosive; the process is non-stationary']);

opts.minProb=.2;
opts.title=sprintf('Event interactions (Hawkes model) %s', condition_name);
plot_tpm_circle(params.branchingRatio', eventnames, opts)

end

