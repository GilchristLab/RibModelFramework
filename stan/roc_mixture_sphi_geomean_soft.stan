/* ============================================================================
 * roc_mixture_sphi_geomean_soft.stan
 *
 * Variant of roc_mixture_sphi_geomean.stan in which the geomean(phi)=1
 * constraint is RELAXED to a soft Gaussian prior:
 *
 *   log(geomean(phi)) = p * mu1 + (1 - p) * mu2
 *   p * mu1 + (1 - p) * mu2 ~ normal(0, geomean_constraint_sd)
 *
 * Both mu1 and mu2 are free parameters here (cf. the hard variant where
 * mu2 is derived).  This adds one more dimension to the posterior, plus
 * a stiff direction in the (mu1, mu2) subspace whose stiffness is
 * inversely proportional to geomean_constraint_sd.
 *
 * Use case: probe how much the hard constraint is fighting the data.
 * Tight sigma (~0.005) -> approximates hard constraint.  Loose sigma
 * (~0.05) -> meaningful give, may reveal a posterior preference for a
 * different log-scale anchor that the hard variant suppresses.
 *
 * Label-switching guard (mu2 >= mu1) uses the same smooth quadratic
 * penalty as the hard variant for consistency.
 *
 * NEW DATA FIELD:
 *   real<lower=0> geomean_constraint_sd     soft prior SD on log(geomean(phi))
 *
 * Default in fit.stan.R: 0.05.
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
    real         mu2_prior_mean;
    real<lower=0> mu2_prior_sd;
    real<lower=0> sigma1_prior_scale;
    real<lower=0> sigma2_prior_scale;

    real<lower=0> geomean_constraint_sd;   // NEW: soft geomean(phi)=1 prior SD

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
    real mu2;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

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

    // Per-gene mixture-LN log prior (threaded)
    target += reduce_sum(partial_mix_prior, gene_indices, grainsize,
                         log_phi, p, mu1, mu2, sigma1, sigma2);

    // Per-gene likelihood (threaded)
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
