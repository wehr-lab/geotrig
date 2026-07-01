%% ── Example script ───────────────────────────────────────────────────────────

T =results(3).TPM;

% --- Euclidean metric ---
[H, H_std, H_vals] = hopkins_statistic(T, 3, 500, 'euclidean');
fprintf('Hopkins H (Euclidean) = %.3f ± %.3f\n', H, H_std);

% --- Symmetric KL divergence (more principled for probability rows) ---
[H_kl, H_kl_std] = hopkins_statistic(T, 3, 500, 'kl_sym');
fprintf('Hopkins H (KL sym)   = %.3f ± %.3f\n', H_kl, H_kl_std);

% --- Interpret ---
interpret_hopkins(H);

% --- Plot distribution of H across trials ---
figure;
histogram(H_vals, 20, 'Normalization', 'probability', ...
          'FaceColor', [0.114 0.620 0.459], 'EdgeColor', 'white');
xline(H, 'r-', sprintf('H = %.3f', H), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(0.5, 'k--', 'random', 'LineWidth', 1);
xlim([0 1]);
xlabel('Hopkins statistic H');
ylabel('Proportion of trials');
title('Hopkins statistic — distribution across Monte-Carlo trials');
grid on; box off;


function interpret_hopkins(H)
    if H < 0.5
        label = 'more uniform than random — no cluster tendency';
    elseif H < 0.75
        label = 'weak cluster tendency';
    elseif H < 0.90
        label = 'moderate cluster tendency';
    else
        label = 'strong cluster tendency';
    end
    fprintf('Interpretation: H = %.3f → %s\n', H, label);
end

% 
% Notes on the implementation Two distance metrics are provided:
% 
% 'euclidean' — treats each row as a point in R6\mathbb{R}^6 R6. Fast and
% standard. 'kl_sym' — symmetric KL divergence
% DKL(p∥q)+DKL(q∥p)D_{KL}(p\|q) + D_{KL}(q\|p) DKL​(p∥q)+DKL​(q∥p). More
% principled for TPM rows since they are probability distributions, not
% arbitrary vectors. Prefer this for your use case.
% 
% Small-nn n considerations: With only 6 states you can't do much better
% than m=2m=2 m=2 or m=3m=3 m=3. The variance will be substantial — hence
% the 500 repetitions. Trust the mean, treat the std as a honesty interval.
% What to expect given your earlier results: The MFPT ratio of 3.55 and
% weak spectral gap suggest HH H will land in the 0.60–0.75 range — real
% but soft cluster tendency — consistent with the picture of leaky,
% non-metastable clusters.
% 
% Thresholds: 
% H<0.5 = more uniform than random; 
% H≈0.5 = no clustering; 
% H>0.75 = moderate tendency;
% H>0.9 = strong. 

% Given your other evidence (ratio 3.55, no spectral gap), you'd expect HH
% H somewhere in the 0.6–0.75 range — real but soft structure.

% output:
% Hopkins H (Euclidean) = 0.454 ± 0.049
% Hopkins H (KL sym)   = 0.299 ± 0.065
% Interpretation: H = 0.454 → more uniform than random — no cluster tendency










