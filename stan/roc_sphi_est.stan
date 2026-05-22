/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC model with sphi estimated.
 *
 * Identical to roc_basic.stan except sphi moves from `data` to `parameters`
 * and gets a half-normal prior (proper, weak).  The lognormal phi prior
 * still uses mphi = -sphi^2/2 (mean(phi) = 1, matching v.3 mphi convention).
 *
 * Prior on sphi: half_normal(0, 5).  Weak and proper.  TODO per user:
 * switch to a flat (improper or uniform on bounded range) prior; uniform
 * on (0, 10) is the easy proper variant.  Mixture-of-lognormals phi prior
 * (multiple sphi components) is a separate Stan file.
 *
 * See roc_basic.stan header for the data layout.
 * ============================================================================ */

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
    real<lower=0> sphi_prior_sd;        // half-normal scale for sphi prior
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;
    real<lower=0> sphi;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // Priors
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi ~ normal(0, sphi_prior_sd);                  // half-normal via lower=0
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);

    // Likelihood (identical to roc_basic.stan)
    for (g in 1:G) {
        target += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);
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
            target += -N_ga[g, a] * log_sum_exp(eta_full);
        }
    }
}
