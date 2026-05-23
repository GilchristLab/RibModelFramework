/* ============================================================================
 * roc_mixture_sphi_ordered.stan
 *
 * Reparameterized variant of roc_mixture_sphi_centered.stan that ELIMINATES
 * both hard rejects from the original by using two Stan idioms:
 *
 *   1. Label-switching guard via Stan's `ordered[2]` type.  We sample
 *      log_component_mean = ordered[2] -- the log of each lognormal
 *      component's MEAN on the natural scale, i.e. log E[phi_k] for k=1,2.
 *      Ordering log_component_mean[1] < log_component_mean[2] guarantees
 *      component 1 is the low-phi (bulk) mode and component 2 is the
 *      high-phi (tail) mode.  Stan's ordered transform is smooth and HMC
 *      handles it natively; no reject() call.
 *
 *   2. mean(phi)=1 constraint as a SOFT prior instead of derived-and-rejected.
 *      Add the constraint log(overall_mean_phi) ~ normal(0, sigma) with
 *      sigma = mean_phi_constraint_sd, supplied via data.  Tight sigma
 *      (e.g. 0.005) approximates the hard constraint to ~0.5%; loose
 *      sigma (e.g. 0.05) gives the sampler more breathing room.
 *
 * The parameter mapping is:
 *   log_component_mean[k] = mu_k + 0.5 * sigma_k^2     (log mean on natural)
 *   mu_k                  = log_component_mean[k] - 0.5 * sigma_k^2
 *
 * The natural-scale mean of the mixture is
 *   E[phi]  = p * exp(log_component_mean[1]) + (1-p) * exp(log_component_mean[2])
 *
 * which we constrain softly toward 1.
 *
 * Semantic match to centered model:
 *   - p (Beta prior) is the prior weight on COMPONENT 1 (low-phi bulk).
 *     Matches native MCMC convention and the existing centered model.
 *   - log_component_mean takes mu1_prior_mean / mu1_prior_sd as its
 *     elementwise prior.  (Same prior on both components; the ordered
 *     transform's Jacobian handles the ordering constraint.)
 *   - sigma1, sigma2: half-normal, scales unchanged.
 *
 * Threading (reduce_sum) and likelihood block are identical to centered.
 *
 * NEW DATA FIELD:
 *   real<lower=0> mean_phi_constraint_sd    soft prior SD on log(E[phi])
 *
 * Default in fit.stan.R: 0.05 (5% deviation tolerated).  For native-
 * equivalent strictness use 0.005 or smaller.
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

    real<lower=0> mean_phi_constraint_sd;   // NEW: soft mean(phi)=1 prior SD

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;
    real<lower=0, upper=1> p;
    ordered[2] log_component_mean;          // log E[phi_k]; built-in ordering
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

    real mu1 = log_component_mean[1] - 0.5 * sigma1 * sigma1;
    real mu2 = log_component_mean[2] - 0.5 * sigma2 * sigma2;
    real overall_mean_phi = p       * exp(log_component_mean[1])
                          + (1 - p) * exp(log_component_mean[2]);
}

model {
    // Hyperpriors
    p                  ~ beta(p_alpha, p_beta);
    log_component_mean ~ normal(mu1_prior_mean, mu1_prior_sd);   // both components
    sigma1             ~ normal(0, sigma1_prior_scale);          // half-normal
    sigma2             ~ normal(0, sigma2_prior_scale);

    // Soft mean(phi)=1 constraint (replaces hard derive-and-reject).
    log(overall_mean_phi) ~ normal(0, mean_phi_constraint_sd);

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
