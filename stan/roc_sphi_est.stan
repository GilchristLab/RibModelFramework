/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC, sphi estimated, centered + non-centered,
 *                      reduce_sum threaded.
 *
 * Unified file replacing roc_sphi_est_centered.stan and
 * roc_sphi_est_noncentered.stan.  The data field `noncentered` selects
 * the parameterization at run time without recompilation:
 *
 *   noncentered = 0  (centered, default):
 *     latent_phi[g] IS log_phi[g].
 *     Prior: log_phi ~ Normal(mphi, sphi)  where mphi = -sphi^2/2
 *     Best when data strongly anchors each gene (G >= 1000, full-genome fits).
 *
 *   noncentered = 1:
 *     latent_phi[g] is z_phi[g] ~ N(0,1).
 *     log_phi[g] = mphi + sphi * z_phi[g]   (transformed parameter)
 *     Removes Neal's funnel in (log_phi, sphi) for data-sparse fits.
 *     Use when centered has poor sphi ESS or chains stuck at small sphi.
 *
 * mphi convention: mphi = -sphi^2/2 (mean(phi)=1, matching v.3 pipeline).
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
 * threads_per_chain > 1 to $sample() to enable reduce_sum parallelism.
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

    int<lower=0, upper=1> noncentered;  // 0=centered (default), 1=non-centered
    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    // latent_phi[g] = log_phi[g] when centered, z_phi[g] when non-centered.
    vector[G] latent_phi;
    real<lower=0> sphi;
}

transformed parameters {
    real mphi = -0.5 * sphi * sphi;
    // log_phi: identity when centered; derived from z_phi when non-centered.
    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi ~ normal(0, sphi_prior_sd);    // half-normal via lower=0 constraint

    // Phi prior: centered form uses the lognormal prior directly on log_phi;
    // non-centered form places std_normal on the latent z_phi.
    if (noncentered) {
        latent_phi ~ std_normal();
    } else {
        latent_phi ~ normal(mphi, sphi);
    }

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
