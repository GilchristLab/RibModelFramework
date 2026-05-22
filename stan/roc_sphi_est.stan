/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC with sphi estimated, non-centered + reduce_sum.
 *
 * Combines:
 *   - sphi as a parameter with half-normal prior (vs roc_basic's data sphi)
 *   - non-centered log_phi parameterization to break Neal's funnel:
 *       z_phi[g] ~ std_normal()                (parameter)
 *       log_phi[g] = -0.5 * sphi^2 + sphi * z_phi[g]   (transformed parameter)
 *   - reduce_sum threaded gene loop for within-chain parallelism
 *
 * Without non-centering, small sphi forces log_phi tightly around -sphi^2/2
 * and HMC step size has to track sphi; sphi ESS_bulk was 91/4000 on the
 * G=1000 centered baseline.
 *
 * THREADING: see roc_basic.stan header.  Compile with STAN_THREADS=true and
 * pass threads_per_chain > 1 to mod$sample() to enable.
 *
 * sphi prior: half_normal(0, sphi_prior_sd).  Switch to flat or uniform via
 * a separate file if needed.  Mixture-of-lognormals (multiple sphi components)
 * is in roc_mixture_sphi.stan.
 *
 * mphi convention: mphi = -sphi^2/2 (mean(phi)=1, matching v.3 convention).
 *
 * See roc_basic.stan header for the data layout.
 * ============================================================================ */

functions {
    real partial_sum_lpdf(array[] int slice_g, int start, int end,
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
    real<lower=0> sphi_prior_sd;        // half-normal scale for sphi prior

    int<lower=1> grainsize;             // reduce_sum partition size
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] z_phi;                    // non-centered latent: z_phi ~ N(0,1)
    real<lower=0> sphi;
}

transformed parameters {
    vector[G] log_phi = -0.5 * sphi * sphi + sphi * z_phi;
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // Priors
    dM    ~ normal(dM_prior_mean, dM_prior_sd);
    dEta  ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi  ~ normal(0, sphi_prior_sd);                 // half-normal via lower=0
    z_phi ~ std_normal();                             // non-centered latent

    // Per-gene likelihood via reduce_sum
    target += reduce_sum(partial_sum_lpdf, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
