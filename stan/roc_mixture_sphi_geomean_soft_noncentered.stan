/* ============================================================================
 * roc_mixture_sphi_geomean_soft_noncentered.stan
 *
 * Non-centered variant of roc_mixture_sphi_geomean_soft.stan.  log_phi is
 * reparametrized via a z-vector centered on component 1:
 *
 *   log_phi[g] = mu1 + sigma1 * z_phi[g]                   (transformed)
 *
 * The mixture prior is then expressed in z-space with component 1 at the
 * standard normal and component 2 translated/scaled relative to it:
 *
 *   z_phi[g] ~ p * N(0, 1) + (1-p) * N(delta, ratio)
 *      delta = (mu2 - mu1) / sigma1
 *      ratio = sigma2 / sigma1
 *
 * Why a non-centered variant.  The centered geomean_soft samples log_phi
 * directly; for small sigma1 the (log_phi, sigma1) plane exhibits a Neal's
 * funnel that HMC navigates poorly.  Reparameterizing log_phi via z_phi
 * removes the funnel in the (z_phi, sigma1) plane for component-1 (bulk)
 * genes.  Component-2 (high-expression) genes still see sigma2 through the
 * `ratio` term, so the funnel in the (z_phi, sigma2) plane is partially
 * mitigated but not eliminated.
 *
 * When to use.  Sparse-data / small-N regimes where the data does not
 * strongly anchor each per-gene log_phi value -- the funnel pathology is
 * most pronounced there.  In data-dense regimes (G in the thousands,
 * 1000+ AA positions per gene) the centered variant should still win on
 * mass-matrix adaptation.  Generally, prefer centered as default and try
 * non-centered if R-hat for sigma1 fails to close at the centered
 * default.
 *
 * The data block, hyperpriors, soft geomean constraint, and likelihood
 * structure are identical to the centered geomean_soft variant.
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

    real partial_mix_prior_z(array[] int slice_g, int start, int end,
                                  vector z_phi,
                                  real p, real delta, real ratio) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += log_mix(
                p,
                normal_lpdf(z_phi[g] | 0, 1),
                normal_lpdf(z_phi[g] | delta, ratio));
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
    real         mu2_prior_mean;
    real<lower=0> mu2_prior_sd;
    real<lower=0> sigma1_prior_scale;
    real<lower=0> sigma2_prior_scale;

    real<lower=0> geomean_constraint_sd;

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
    vector[G] z_phi;                       // non-centered latent for log_phi
    real<lower=0, upper=1> p;
    real mu1;
    real mu2;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    // Centered-on-component-1 non-centering.
    vector[G] log_phi = mu1 + sigma1 * z_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    // Component-2 location and scale in z-space.
    real delta = (mu2 - mu1) / sigma1;
    real ratio = sigma2 / sigma1;

    real log_geomean_phi = p * mu1 + (1.0 - p) * mu2;
}

model {
    // Hyperpriors
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    mu2    ~ normal(mu2_prior_mean, mu2_prior_sd);
    sigma1 ~ normal(0, sigma1_prior_scale);
    sigma2 ~ normal(0, sigma2_prior_scale);

    // Soft geomean(phi)=1 constraint
    log_geomean_phi ~ normal(0, geomean_constraint_sd);

    // Label-switching guard (smooth quadratic penalty).
    target += -PENALTY_STRENGTH * square(fmax(0, mu1 - mu2));

    // dM, dEta priors
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Mixture prior is on log_phi in z-space (per-component normal in z).
    // Note: log_phi = mu1 + sigma1 * z_phi, and the prior on z_phi is the
    // shifted/scaled mixture above.  This is the standard non-centered
    // reparametrization where the sampled parameter (z_phi) carries the
    // mixture prior; log_phi is a deterministic function.  No explicit
    // Jacobian is needed because z_phi (not log_phi) is the parameter.
    target += reduce_sum(partial_mix_prior_z, gene_indices, grainsize,
                         z_phi, p, delta, ratio);

    // Per-gene likelihood (threaded; depends on log_phi / phi)
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
