/* ============================================================================
 * roc_mixture_sphi.stan -- ROC HMC with 2-component lognormal phi prior,
 *                          non-centered (component-1 anchored) + reduce_sum.
 *
 * Combines:
 *   - 2-component lognormal mixture prior on phi (mirrors C++
 *     Parameter::deriveMu2 with PHI_CONSTRAINT_MEAN=1)
 *   - non-centered log_phi anchored on the dominant component 1:
 *       z_phi[g] ~ std_normal()                (parameter)
 *       log_phi[g] = mu1 + sigma1 * z_phi[g]   (transformed parameter)
 *     Component 1's contribution to the mixture becomes std_normal_lpdf(z_phi)
 *     exactly; component 2's contribution retains a small residual funnel
 *     (rare in posterior draws since E[p] = 0.8 under Beta(8, 2) prior).
 *   - reduce_sum threaded gene loop for within-chain parallelism
 *
 * Change-of-variables math: with log_phi = mu1 + sigma1 * z_phi, the original
 * mixture-LN prior on log_phi
 *
 *   log_mix(p, N(log_phi | mu1, sigma1), N(log_phi | mu2, sigma2))
 *
 * becomes in z_phi-space (|d log_phi / d z_phi| = sigma1):
 *
 *   log_mix(p,
 *           std_normal_lpdf(z_phi),
 *           N(log_phi | mu2, sigma2) + log(sigma1))
 *
 * The log(sigma1) Jacobian appears ONLY on the component-2 leg because it
 * cancels exactly on the component-1 leg.
 *
 * mu2 DERIVED from PHI_CONSTRAINT_MEAN (= mean(phi) = 1):
 *
 *   p * exp(mu1 + sigma1^2/2) + (1 - p) * exp(mu2 + sigma2^2/2) = 1
 *   => mu2 = log( (1 - p * exp(mu1 + sigma1^2/2)) / (1 - p) ) - sigma2^2/2
 *
 * Label-switching guard: mu2 >= mu1 (component 2 has the higher mode in
 * log-space).  Constraint infeasibility (numer <= 0) -> reject.
 *
 * THREADING: see roc_basic.stan header.  Compile with STAN_THREADS=true and
 * pass threads_per_chain > 1 to mod$sample() to enable.
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

    // Per-gene mixture-LN prior contribution in z_phi-space.
    // Component 1 leg -> std_normal_lpdf(z_phi[g]);
    // component 2 leg -> normal_lpdf(log_phi[g] | mu2, sigma2) + log(sigma1).
    // Combined via log_mix(p, leg1, leg2).
    real partial_mix_prior_lpdf(array[] int slice_g, int start, int end,
                                vector z_phi, vector log_phi,
                                real p, real mu2, real sigma1, real sigma2) {
        real lp = 0;
        real log_sigma1 = log(sigma1);
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += log_mix(
                p,
                std_normal_lpdf(z_phi[g]),
                normal_lpdf(log_phi[g] | mu2, sigma2) + log_sigma1);
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
    vector[G] z_phi;                        // non-centered latent for component 1
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector[G] log_phi = mu1 + sigma1 * z_phi;
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
    target += reduce_sum(partial_mix_prior_lpdf, gene_indices, grainsize,
                         z_phi, log_phi, p, mu2, sigma1, sigma2);

    // Per-gene likelihood (threaded)
    target += reduce_sum(partial_sum_lpdf, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
