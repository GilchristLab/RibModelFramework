/* ============================================================================
 * roc_mixture_sphi_centered_smooth.stan
 *
 * Variant of roc_mixture_sphi_centered.stan that replaces the hard reject()
 * calls (mu2<mu1 label-switching guard; numer<=0 mean-constraint
 * infeasibility) with SMOOTH quadratic penalties.  Hard rejects break HMC:
 * any leapfrog trajectory that crosses the boundary returns -inf log-prob,
 * the whole NUTS subtree is discarded, mass-matrix adaptation is corrupted,
 * and divergences inflate.  Quadratic penalties have gradient zero at the
 * boundary (smooth) but grow steeply outside, pushing trajectories back.
 *
 * The constraints themselves are unchanged in spirit:
 *   numer = 1 - p * exp(mu1 + sigma1^2/2)  must be > 0 for mu2 to be defined
 *   mu2 >= mu1 to prevent label-switching
 *
 * Penalty strength PENALTY_STRENGTH = 1e3 was chosen so that boundary
 * violations of magnitude ~0.1 incur penalties of ~10 log-prob units,
 * comparable to typical likelihood-scale gradients on this model.
 * Adjust if you see boundary-bleed or, conversely, persistent cliff-like
 * behavior on the boundary.
 *
 * mu2 derivation uses fmax(raw_numer, NUMER_EPS) to keep log() defined
 * even when raw_numer briefly goes negative during a trajectory (the
 * penalty term will pull the chain back).
 *
 * Everything else (data block, likelihood, hyperpriors, reduce_sum) is
 * identical to roc_mixture_sphi_centered.stan.
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

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;

    // Smooth-penalty tuning constants.  See header for rationale.
    real PENALTY_STRENGTH = 1.0e3;
    real NUMER_EPS        = 1.0e-8;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

    // raw_numer can briefly go negative during a leapfrog trajectory;
    // we use safe_numer (clipped to >= NUMER_EPS) for log() definition,
    // while raw_numer drives the smooth penalty so the chain is pulled
    // back to the feasible region.
    real raw_numer  = 1.0 - p * exp(mu1 + 0.5 * sigma1 * sigma1);
    real safe_numer = fmax(raw_numer, NUMER_EPS);
    real mu2 = log(safe_numer / (1.0 - p)) - 0.5 * sigma2 * sigma2;
}

model {
    // Hyperpriors (unchanged from centered model)
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    sigma1 ~ normal(0, sigma1_prior_scale);
    sigma2 ~ normal(0, sigma2_prior_scale);

    // SMOOTH PENALTIES (replace hard reject() from centered model).
    // square(fmax(0, x)) has gradient zero at x=0 (smooth from both sides),
    // grows quadratically for x > 0 -- ideal for HMC.
    target += -PENALTY_STRENGTH * square(fmax(0, -raw_numer));    // raw_numer >= 0
    target += -PENALTY_STRENGTH * square(fmax(0, mu1 - mu2));     // mu2 >= mu1

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
