/* ============================================================================
 * roc_mixture_sphi_geomean.stan
 *
 * Variant of the mixture-LN ROC model that anchors the log-scale via a
 * GEOMETRIC-MEAN constraint instead of the arithmetic-mean constraint
 * used by roc_mixture_sphi_centered{,_smooth}.stan.
 *
 * The geomean(phi) = 1 constraint is:
 *   log(geomean(phi)) = E[log phi] = p * mu1 + (1 - p) * mu2 = 0
 *
 * This is LINEAR in (mu1, mu2), which means:
 *   - It is always satisfiable (no "numer <= 0" boundary).
 *   - mu2 is derived from mu1 as a closed-form linear function:
 *         mu2 = -p * mu1 / (1 - p)
 *   - The constraint manifold is a flat hyperplane in the (mu1, mu2)
 *     subspace, well-conditioned for HMC.
 *
 * Compared to the arithmetic-mean (centered_smooth) variant:
 *   - One fewer free hyperparameter (mu2 is fully determined by mu1, p).
 *   - No safe_numer clipping, no NUMER_EPS, no "infeasibility" penalty --
 *     only the label-switching guard remains.
 *   - Phi values are anchored differently: median(phi) is closer to 1 than
 *     under mean=1 (where mean(phi)=1 forces median to drift below 1 for
 *     right-skewed mixtures).
 *
 * Label-switching guard (mu2 >= mu1) is enforced via the same smooth
 * quadratic penalty as the smooth variant.  Hard rejects are avoided.
 *
 * Threading (reduce_sum) and likelihood block are identical to centered.
 *
 * NOTE on biological equivalence:
 *   A fit under geomean=1 and a fit under mean=1 describe the SAME
 *   biological process; they differ only in normalization of the log
 *   phi scale.  Under geomean=1 the natural-scale mean(phi) > 1 (since
 *   E[exp(log_phi)] >= exp(E[log_phi]) by Jensen, with equality only when
 *   the distribution is degenerate).  Specifically:
 *       mean(phi) under geomean=1 = exp((p*sigma1^2 + (1-p)*sigma2^2)/2
 *                                       + (mixing variance term))
 *
 * ============================================================================ */

functions {
    real partial_sum(array[] int slice_g, int start, int end,
                          int A,
                          array[] int aa_start, array[] int aa_end,
                          array[,] int y_k,
                          array[,] int N_ga,
                          vector dM, vector dEta,
                          vector phi) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);
            for (a in 1:A) {
                if (N_ga[g, a] == 0) continue;
                int s = aa_start[a];
                int e = aa_end[a];
                int n = e - s + 1;
                vector[n + 1] eta_full;
                eta_full[1] = 0;
                for (k in 1:n) {
                    eta_full[k + 1] = -dM[s - 1 + k] - dEta[s - 1 + k] * phi[g];
                }
                lp += -N_ga[g, a] * log_sum_exp(eta_full);
            }
        }
        return lp;
    }

    real partial_mix_prior(array[] int slice_g, int start, int end,
                                vector log_phi,
                                real p, real mu1, real mu2,
                                real sigma1, real sigma2) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += log_mix(
                p,
                normal_lpdf(log_phi[g] | mu1, sigma1),
                normal_lpdf(log_phi[g] | mu2, sigma2));
        }
        return lp;
    }
}

data {
    int<lower=1> G;
    int<lower=1> A;
    int<lower=1> K;
    array[A] int<lower=1> aa_start;
    array[A] int<lower=1> aa_end;
    array[G, K] int<lower=0> y_k;
    array[G, A] int<lower=0> N_ga;
    vector[K] dM_prior_mean;
    vector<lower=0>[K] dM_prior_sd;
    vector[K] dEta_prior_mean;
    vector<lower=0>[K] dEta_prior_sd;

    real<lower=0> p_alpha;
    real<lower=0> p_beta;
    real         mu1_prior_mean;
    real<lower=0> mu1_prior_sd;
    real<lower=0> sigma1_prior_scale;
    real<lower=0> sigma2_prior_scale;

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;

    real PENALTY_STRENGTH = 1.0e3;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

    // Geomean=1 hard linear constraint: p*mu1 + (1-p)*mu2 = 0
    // => mu2 = -p * mu1 / (1 - p)
    // Always defined (1-p > 0 strictly via Beta prior support).
    real mu2 = -p * mu1 / (1.0 - p);
}

model {
    // Hyperpriors
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    sigma1 ~ normal(0, sigma1_prior_scale);
    sigma2 ~ normal(0, sigma2_prior_scale);

    // Label-switching guard (smooth quadratic penalty; same as smooth variant).
    // mu1 should be below mu2 for component-1 = bulk-low, component-2 = tail-high.
    target += -PENALTY_STRENGTH * square(fmax(0, mu1 - mu2));

    // dM, dEta priors
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Per-gene mixture-LN log prior (threaded)
    target += reduce_sum(partial_mix_prior, gene_indices, grainsize,
                         log_phi, p, mu1, mu2, sigma1, sigma2);

    // Per-gene likelihood (threaded)
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
