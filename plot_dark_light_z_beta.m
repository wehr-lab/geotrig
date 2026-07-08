% plot dark/light beta for engagement drift for each bin

%run analyze_localized_laser_effects first to generate stuff

laserBeta = laserResults.beta;      % effect of laser-on, per bin, in dark
laserSE   = laserResults.SE;
lightBeta = -darkResults.beta;      % NEGATE: darkResults.beta is dark-vs-light;
lightSE   = darkResults.SE;         % light-vs-dark is its negation. SE unaffected by sign flip.
binCenters = laserResults.binCenters;   % same grid in both runs by construction

% fprintf('\n  %-6s %-8s %-14s %-14s %-14s\n', 'bin', 'z', 'laser beta', 'light beta', 'diff (z-test)');
% diffZ = nan(size(binCenters)); diffP = nan(size(binCenters));
% for b = 1:numel(binCenters)
%     if isnan(laserBeta(b)) || isnan(lightBeta(b)), continue; end
%     d = laserBeta(b) - lightBeta(b);
%     seD = sqrt(laserSE(b)^2 + lightSE(b)^2);
%     diffZ(b) = d / seD;
%     diffP(b) = 2*(1 - normcdf(abs(diffZ(b))));
%     fprintf('  %-6d %-8.2f %-14.4f %-14.4f z=%.2f, p=%.4g\n', ...
%         b, binCenters(b), laserBeta(b), lightBeta(b), diffZ(b), diffP(b));
% end

% validBoth = ~isnan(laserBeta) & ~isnan(lightBeta);
% [rho, pRho] = corr(laserBeta(validBoth), lightBeta(validBoth), 'Type', 'Pearson');
% fprintf('\n  Cross-bin correlation (laser beta vs light beta): r = %.3f, p = %.4g (n = %d bins)\n', ...
%     rho, pRho, nnz(validBoth));
% fprintf('  If the mirror hypothesis holds: expect r close to +1 (same shape, same sign,\n');
% fprintf('  similar magnitude) and few/no bins with a significant diff (z-test) above.\n');
% fprintf('  A weak/non-significant r, or a cluster of significant per-bin differences,\n');
% fprintf('  is evidence the two manipulations act on different parts of the range.\n');

figure('Name', 'light effect, per bin');
plot(binCenters, lightBeta, '-o', 'DisplayName', 'light-vs-dark effect (laser off)');
yline(0, '--k', 'HandleVisibility','off');
xlabel('z (engagement axis)', 'FontSize', 18); ylabel('per-bin drift beta', 'FontSize', 18);
title('light-vs-dark per-bin effect ', 'FontSize', 18);
ylim([-1.1 1.6])

    local_overlayStateRanges(z, state_id, inSubset);
