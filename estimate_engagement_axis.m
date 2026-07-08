%% =========================================================================
%  ESTIMATE_ENGAGEMENT_AXIS
%  Bottom-up estimation of a continuous engagement-intensity / range axis
%  from kinematic features, via PCA on sign-oriented, z-scored features.
%
%  Features: mouse_spd, range, az (azimuth), cricket_spd
%  Output:   z(t) = PC1, oriented so that higher z = higher engagement
%            (fast, close range, on-target heading)
% =========================================================================

%% 0. ASSUMPTIONS -- check these against your actual data before trusting z(t)
% -------------------------------------------------------------------------
% mouse_spd   : Nx1, predator speed, higher = faster
% range       : Nx1, predator-prey distance, LOWER = closer/more engaged
% az          : Nx1, azimuth of prey relative to predator heading
%               ASSUMED CONVENTION: az = 0 means prey dead-ahead;
%               larger |az| = prey more off-axis (poorer target lock).
%               If your az is instead a signed lateral bearing where the
%               *sign* carries information you want to keep (e.g. left vs
%               right turning bias), do NOT collapse to |az| -- use az
%               directly and let PCA/sign-correction handle it, or split
%               it into a separate "turning" axis rather than folding it
%               into engagement.
% cricket_spd : Nx1, prey speed
%
% All vectors must be the same length, frame-aligned, with NaN wherever
% data is missing or the frame is NOTA/censored.

az_unit = 'deg';   % 'deg' or 'rad' -- set to match your az variable

%% 1. ASSEMBLE FEATURES, RESOLVE AZIMUTH CIRCULARITY -----------------------
if strcmp(az_unit, 'deg')
    az_rad = deg2rad(az);
else
    az_rad = az;
end

% Wrap to (-pi, pi], take magnitude -> 0 = on-target, pi = facing away.
% (Mapping Toolbox not required -- this is a manual wrap.)
az_dev = abs(mod(az_rad + pi, 2*pi) - pi);

X_raw = [mouse_spd(:), range(:), az_dev(:), cricket_spd(:)];
feat_names = {'mouse_spd', 'range', 'az_dev', 'cricket_spd'};

%% 2. NaN-AWARE Z-SCORING ---------------------------------------------------
% NOTE: if you have multiple animals/sessions with different baseline
% speeds or arena sizes, consider z-scoring WITHIN group instead of
% globally -- see grouped-normalization snippet at the bottom of this file.
mu_x = mean(X_raw, 1, 'omitnan');
sd_x = std(X_raw, 0, 1, 'omitnan');
Xz   = (X_raw - mu_x) ./ sd_x;

%% 3. ORIENT FEATURES SO HIGHER VALUE = HIGHER ENGAGEMENT -------------------
% range and az_dev both DECREASE with engagement -> flip sign
sign_flip = [+1, -1, -1, +1];   % [mouse_spd, range, az_dev, cricket_spd]
Xo = Xz .* sign_flip;

%% 4. PCA ON COMPLETE ROWS (drop NaN / NOTA frames) -------------------------
valid = all(~isnan(Xo), 2);
if nnz(valid) < 100
    warning('Fewer than 100 complete-feature frames (%d) -- check upstream NaNs.', nnz(valid));
end

[coeff, score, ~, ~, explained] = pca(Xo(valid, :));

fprintf('Variance explained by PC1: %.1f%%\n', explained(1));
fprintf('PC1 loadings (sign-oriented features):\n');
for i = 1:numel(feat_names)
    fprintf('  %-14s % .3f\n', feat_names{i}, coeff(i, 1));
end
% All loadings should be positive and roughly comparable in magnitude if
% the sign-orientation step worked and no single feature is dominating
% for a boring reason (e.g. a near-constant, low-variance channel).
% A loading near zero for one feature suggests it isn't tracking this
% axis and may deserve its own dimension (see Section 8).

%% 5. SIGN-CORRECT PC1 AGAINST AN UNAMBIGUOUS ANCHOR -------------------------
% PCA sign is arbitrary; anchor against oriented range (closer = higher
% engagement is the least contestable of the four features).
anchor = Xo(valid, 2);
if corr(score(:, 1), anchor, 'rows', 'complete') < 0
    score(:, 1) = -score(:, 1);
    coeff(:, 1) = -coeff(:, 1);
end

%% 6. REINSERT INTO FULL-LENGTH VECTOR (NaN where excluded) -----------------
z = nan(size(range));
z(valid) = score(:, 1);

%% 7. OPTIONAL SMOOTHING -----------------------------------------------------
win = 5;  % frames -- tune to your frame rate (200 fps) and feature noise
z_smooth = movmean(z, win, 'omitnan');

%% 8. VALIDATION AGAINST DISCRETE STATE LABELS -------------------------------
% Requires a parallel Nx1 `state_id` coded 1-6 in the order:
% statenames = {'hot pursuit','chase','following','stalk','wander','pause'};

% ── Build stateMask% ──────────────────────────────────────────────(mike)
stateMask = logical([hotpursuit(:), chase(:), follow(:), stalk(:), wander(:), pause(:)]);
nStates=size(stateMask, 2);
num_frames = size(stateMask, 1);
stateSeq=zeros(num_frames,1);
for s = 1:nStates
    stateSeq(stateMask(:,s)) = s;
end
state_id=stateSeq;

if exist('state_id', 'var')
    statenames = {'hot pursuit', 'chase', 'following', 'stalk', 'wander', 'pause'};


    % Exclude NOTA (state_id == 0) frames -- boxplot/findgroups need a
    % clean 1:6 categorical, and NOTA isn't part of the 6-state ordering.
    keep = state_id >= 1 & state_id <= 6 & ~isnan(z);

    figure;
    boxplot(z(keep), state_id(keep), 'Labels', statenames);
    ylabel('Latent engagement axis, z (PC1)');
    title('Kinematic latent axis by behavioral state');
    grid on;

    g = findgroups(state_id(keep));
    med_by_state = splitapply(@(x) median(x, 'omitnan'), z(keep), g);

    fprintf('\nMedian z by state (expect monotonic decrease, state 1 -> 6):\n');
    for s = 1:numel(statenames)
        fprintf('  %-12s % .3f\n', statenames{s}, med_by_state(s));
    end

    if all(diff(med_by_state) < 0)
        disp('Monotonic ordering confirmed: PCA axis agrees with the assumed state ordering.');
    else
        warning(['Monotonic ordering NOT confirmed. Check: (a) az sign convention, ' ...
            '(b) whether one feature is dominating the loadings, or ' ...
            '(c) whether the discrete state ordering itself needs revisiting.']);
    end
end
 

%% 9. PLOT EXAMPLE TRACE ------------------------------------------------------
figure;
plot(z_smooth, 'LineWidth', 1);
xlabel('Frame'); ylabel('z(t) -- engagement intensity axis');
title('Continuous latent engagement trajectory');
grid on;

%% =========================================================================
%  OPTIONAL: grouped (per-animal/session) z-scoring instead of Section 2
%  Use this if baseline speeds/arena scale differ meaningfully across
%  animals or sessions. Requires a parallel Nx1 `group_id` vector.
% =========================================================================
% g = findgroups(group_id);
% Xz = zeros(size(X_raw));
% for k = 1:max(g)
%     idx = (g == k);
%     mu_k = mean(X_raw(idx,:), 1, 'omitnan');
%     sd_k = std(X_raw(idx,:), 0, 1, 'omitnan');
%     Xz(idx,:) = (X_raw(idx,:) - mu_k) ./ sd_k;
% end
