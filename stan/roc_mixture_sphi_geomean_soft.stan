/* ============================================================================
 * roc_mixture_sphi_geomean_soft.stan
 *
 * Unified mixture-LN ROC model with soft geomean(phi)=1 constraint and
 * optional non-centered reparametrization.  The data field `noncentered`
 * selects the parameterization without recompilation:
 *
 *   noncentered = 0  (centered, default):
 *     latent_phi[g] IS log_phi[g].
 *     Mixture prior: log_phi ~ p*N(mu1,sigma1) + (1-p)*N(mu1+sep,sigma2)
 *     Best when data strongly anchors each gene (G >= 1000).
 *
 *   noncentered = 1:
 *     latent_phi[g] is z_phi[g] anchored on component 1.
 *     log_phi[g] = mu1 + sigma1 * z_phi[g]   (transformed parameter)
 *     Mixture prior in z-space: z_phi ~ p*N(0,1) + (1-p)*N(delta, ratio)
 *       where delta = sep/sigma1, ratio = sigma2/sigma1
 *     Removes (log_phi, sigma1) Neal's funnel for the bulk component.
 *     Use when centered has E-BFMI < 0.3 or R-hat > 1.1 on sigma1.
 *
 * Label-switching prevention: mu2 is reparametrized as mu1 + sep where
 * sep > 0 is a constrained parameter with an informative prior.  This
 * replaces the old quadratic penalty (which had zero gradient in the valid
 * region and failed to prevent warmup mode-flips).
 *
 * Soft geomean constraint in both variants:
 *   log(geomean(phi)) = p*mu1 + (1-p)*mu2 ~ N(0, geomean_constraint_sd)
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
 * threads_per_chain > 1 to enable reduce_sum parallelism.
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

    // Centered mixture prior on log_phi.
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

    // Non-centered mixture prior on z_phi (component 1 anchored at N(0,1)).
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
    real         sep_prior_mean;
    real<lower=0> sep_prior_sd;
    real<lower=0> sigma1_prior_scale;
    real<lower=0> sigma2_prior_scale;

    real<lower=0> geomean_constraint_sd;

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
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sep;   // mu2 - mu1; constrained positive to prevent label-switching
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    real mu2 = mu1 + sep;

    // log_phi: identity when centered; derived from z_phi when non-centered.
    vector[G] log_phi = noncentered ? (mu1 + sigma1 * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    real log_geomean_phi = mu1 + (1.0 - p) * sep;

    // Component-2 location and scale in z-space (used only when noncentered=1;
    // computed unconditionally to avoid branching in transformed parameters).
    real delta = sep / sigma1;
    real ratio = sigma2 / sigma1;
}

model {
    // Hyperpriors
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    sep    ~ normal(sep_prior_mean, sep_prior_sd);  // <lower=0> enforces mu2>mu1
    sigma1 ~ normal(0, sigma1_prior_scale);
    sigma2 ~ normal(0, sigma2_prior_scale);

    // Soft geomean(phi)=1 constraint
    log_geomean_phi ~ normal(0, geomean_constraint_sd);

    // dM, dEta priors
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Mixture prior: centered form uses log_phi directly; non-centered uses
    // z_phi with component-2 shifted/scaled by (delta, ratio).
    if (noncentered) {
        target += reduce_sum(partial_mix_prior_z, gene_indices, grainsize,
                             latent_phi, p, delta, ratio);
    } else {
        target += reduce_sum(partial_mix_prior, gene_indices, grainsize,
                             latent_phi, p, mu1, mu2, sigma1, sigma2);
    }

    // Per-gene likelihood (threaded)
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
