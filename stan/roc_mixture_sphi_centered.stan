/* ============================================================================
 * roc_mixture_sphi_centered.stan -- ROC HMC with 2-component lognormal phi
 *                                   prior, CENTERED + reduce_sum.
 *
 * log_phi is a direct parameter with the 2-component lognormal mixture prior
 * evaluated directly via log_mix.  No reparameterization on log_phi.  This
 * is the data-informative workhorse companion to roc_sphi_est_centered.stan;
 * see that file's header for why centered is the default at typical G.
 *
 * Mirrors the phi-mixture-LN prior implemented in C++ on
 * feat/rst-mixture-state (task #12 family).  mu2 DERIVED from
 * PHI_CONSTRAINT_MEAN = mean(phi) = 1:
 *
 *   p * exp(mu1 + sigma1^2/2) + (1 - p) * exp(mu2 + sigma2^2/2) = 1
 *   => mu2 = log( (1 - p * exp(mu1 + sigma1^2/2)) / (1 - p) ) - sigma2^2/2
 *
 * Label-switching guard: mu2 >= mu1 (component 2 is the higher-mode lognormal
 * in log-space).  Constraint infeasibility (numer <= 0) -> reject.
 *
 * THREADING: both the per-gene mixture log-prior and the per-gene likelihood
 * are partitioned via reduce_sum.  Compile with cpp_options=list(stan_threads
 * =TRUE) and pass threads_per_chain > 1 to enable.
 *
 * Hyperpriors match C++ defaults in src/Parameter.cpp:
 *   p      ~ Beta(p_alpha, p_beta)         default Beta(8, 2)
 *   mu1    ~ Normal(mu1_prior_mean, mu1_prior_sd)   default N(0, 10)
 *   sigma1 ~ half_normal(0, sigma1_prior_scale)     default scale 1
 *   sigma2 ~ half_normal(0, sigma2_prior_scale)     default scale 1
 *
 * See roc_basic.stan header for the data layout (identical here; only the
 * prior on log_phi changes).
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

    // Per-gene mixture-LN log prior, evaluated directly at log_phi[g]:
    //   log_mix(p, N(log_phi | mu1, sigma1), N(log_phi | mu2, sigma2))
    // No Jacobian needed (centered form).
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

    int<lower=1> grainsize;             // reduce_sum partition size
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;                  // CENTERED: log_phi is a direct parameter
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

    // Derived mu2 from mean=1 constraint
    real numer = 1.0 - p * exp(mu1 + 0.5 * sigma1 * sigma1);
    real mu2 = log(numer / (1.0 - p)) - 0.5 * sigma2 * sigma2;
}

model {
    // Hyperpriors
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    sigma1 ~ normal(0, sigma1_prior_scale);             // half-normal via lower=0
    sigma2 ~ normal(0, sigma2_prior_scale);

    // Constraint feasibility + label-switching guard
    if (numer <= 0) reject("mixture-LN constraint infeasible: numer <= 0");
    if (mu2 < mu1)  reject("label-switching guard: mu2 < mu1");

    // dM, dEta priors (Gaussian, per-codon means/SDs)
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
