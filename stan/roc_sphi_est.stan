/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC model with sphi estimated, NON-CENTERED on log_phi.
 *
 * Identical to roc_basic.stan except (a) sphi moves from `data` to `parameters`
 * with a half-normal prior (proper, weak), and (b) log_phi is non-centered:
 *
 *   z_phi[g] ~ std_normal()                (parameter)
 *   log_phi[g] = -0.5 * sphi^2 + sphi * z_phi[g]    (transformed parameter)
 *
 * This breaks Neal's funnel between sphi and the per-gene log_phi: with the
 * old centered form (`log_phi ~ normal(-sphi^2/2, sphi)`), small sphi forces
 * log_phi tightly around -sphi^2/2 and HMC step size has to track sphi.
 * Non-centering decouples z_phi from sphi so the geometry is sphi-invariant.
 *
 * Earlier centered version hit sphi ESS_bulk = 91 / 4000 on sim S288c
 * (G=1000, sphi=1.0).  This non-centered form should restore sphi ESS to
 * the per-iteration mixing level of the other scalar parameters.
 *
 * mphi convention: mphi = -sphi^2/2 (mean(phi)=1, matching v.3 convention).
 * sphi prior: half_normal(0, sphi_prior_sd).  Switch to flat or uniform via
 * a separate file if needed.
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
    vector[G] z_phi;                    // non-centered latent: z_phi ~ N(0,1)
    real<lower=0> sphi;
}

transformed parameters {
    // log_phi is deterministic in (sphi, z_phi); funnel-free in the sampling
    // space because z_phi has unit-scale geometry independent of sphi.
    vector[G] log_phi = -0.5 * sphi * sphi + sphi * z_phi;
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // Priors
    dM    ~ normal(dM_prior_mean, dM_prior_sd);
    dEta  ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi  ~ normal(0, sphi_prior_sd);                 // half-normal via lower=0
    z_phi ~ std_normal();                             // non-centered latent

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
