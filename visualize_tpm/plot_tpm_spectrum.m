function eigvals = plot_tpm_spectrum(P, varargin)
% PLOT_TPM_SPECTRUM  Plot eigenvalues of a TPM in the complex plane (with
% the unit circle for reference), as a diagnostic for natural timescale
% separation / number of "slow" coarse states.
%
% USAGE:
%   eigvals = plot_tpm_spectrum(P)
%   eigvals = plot_tpm_spectrum(P, 'nLabel', 5)
%
% INPUTS:
%   P - N x N transition probability matrix (rows sum to 1)
%
% OPTIONAL NAME-VALUE ARGS:
%   'nLabel' - number of largest-magnitude eigenvalues to label with
%              their index and value (default: min(N,5))
%
% OUTPUT:
%   eigvals - N x 1 vector of eigenvalues, sorted by descending magnitude
%
% NOTES:
%   - All eigenvalues lie within the unit circle (|lambda| <= 1), with
%     lambda = 1 always present for an irreducible chain.
%   - A visible GAP between |lambda_1|=1 and the next eigenvalue's
%     magnitude (the "spectral gap") suggests fast mixing overall.
%   - If there's a cluster of eigenvalues near 1 (e.g., |lambda_2|, ...,
%     |lambda_K| all close to 1) separated by a gap from the rest, this
%     suggests K natural "slow" coarse states / metastable clusters --
%     a useful starting point for choosing K in your lumping.
%   - A second plot shows |lambda| sorted in descending order, which
%     often makes such gaps easier to spot than the complex-plane view.

p = inputParser;
addRequired(p, 'P', @(x) ismatrix(x) && size(x,1)==size(x,2));
addParameter(p, 'nLabel', [], @(x) isempty(x) || (isscalar(x) && x>=0));
parse(p, P, varargin{:});

fps=200;

P = p.Results.P;
N = size(P,1);

nLabel = p.Results.nLabel;
if isempty(nLabel)
    nLabel = min(N,5);
end

eigvals = eig(P);
[~, order] = sort(abs(eigvals), 'descend');
eigvals = eigvals(order);

figure('Color','w');

% --- Left panel: complex plane ---
subplot(1,3,1);
hold on;
th = linspace(0,2*pi,200);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1); % unit circle
plot(real(eigvals), imag(eigvals), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k');

for k = 1:min(nLabel, N)
    lam = eigvals(k);
    text(real(lam)+0.03, imag(lam)+0.03, sprintf('\\lambda_{%d}', k), ...
        'FontSize', 10);
end

axis equal;
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title('Eigenvalues of P (complex plane)');
grid on;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
hold off;

% --- middle panel: |lambda| sorted, descending ---
subplot(1,3,2);
mags = abs(eigvals);
bar(mags, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Index k (sorted by |\lambda_k|, descending)');
ylabel('|\lambda_k|');
title('Eigenvalue Magnitudes');
ylim([0 1.05]);
grid on;

for k = 1:min(nLabel, N)
    text(k, mags(k)+0.03, sprintf('%.3f', mags(k)), ...
        'HorizontalAlignment','center', 'FontSize', 9);
end

sgtitle('TPM Spectrum: look for gaps suggesting natural cluster counts');

% --- right panel:     % Relaxation timescales (skip eigenvalue 1 = stationary)

subplot(1,3,3);
mags = abs(eigvals);
bar(mags, 'FaceColor', [0.2 0.4 0.8]);

tauRelax = -1 ./ log(mags(2:end));
tauRelax(isinf(tauRelax)) = NaN;

bar(2:length(tauRelax)+1, tauRelax/fps, 'FaceColor', [0.2 0.4 0.8]);

xlabel('eigenvalue #');
ylabel('Relaxation timescale (sec)');
title('Implied timescales  \tau = -1/log|\lambda|'); grid on;

set(gcf, 'pos', [ 614   818   946   420])
end
