/* ============================================================================
 * roc_sphi_est_centered.stan -- ROC HMC, sphi estimated, CENTERED + reduce_sum.
 *
 * log_phi is a direct parameter with prior log_phi ~ Normal(-sphi^2/2, sphi).
 * sphi is a parameter with half-normal(0, sphi_prior_sd) prior.  No funnel
 * reparameterization.  This is the "data-informative" workhorse: when the
 * data strongly informs each log_phi[g] (which is the typical regime for
 * full-genome fits with G >= 1000 and standard codon counts per gene), the
 * centered parameterization mixes much better than non-centered because
 * the data pins log_phi tightly and the (sphi, log_phi) joint is not
 * dominated by the prior.
 *
 * The non-centered version (roc_sphi_est_noncentered.stan) helps only in
 * the data-sparse regime (very small G, or genes with few codons), where
 * the prior dominates and small sphi creates a centered-form funnel.
 *
 * mphi convention: mphi = -sphi^2/2 (mean(phi)=1, matching v.3 convention).
 *
 * THREADING: per-gene likelihood via reduce_sum.  Compile with
 * cpp_options=list(stan_threads=TRUE) and pass threads_per_chain > 1 to
 * enable.
 *
 * See roc_basic.stan header for the data layout.
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
    real<lower=0> sphi;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // Priors
    dM      ~ normal(dM_prior_mean, dM_prior_sd);
    dEta    ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi    ~ normal(0, sphi_prior_sd);                 // half-normal via lower=0
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);         // lognormal phi prior (mean=1)

    // Per-gene likelihood via reduce_sum
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
