/* ============================================================================
 * roc_basic.stan -- ROC model in Stan, reduce_sum threaded version.
 *
 * Per-AA multinomial-logit on codon counts: for each gene g and AA a,
 *
 *   eta_{g,a,c} = 0                              if c is the reference codon
 *               = -dM_c - dEta_c * phi_g         otherwise
 *
 *   P(codon c | a, g) = exp(eta_{g,a,c}) / sum_{c' in a} exp(eta_{g,a,c'})
 *
 *   y_{g,a,:}  ~  multinomial(P(. | a, g))
 *
 * Parameters jointly sampled by HMC: dM (K non-ref codons), dEta (K),
 * log_phi (G genes).  Reference codons have dM = dEta = 0 by construction
 * (not free parameters).  sphi is held fixed for this version, with the
 * lognormal phi prior constraint mphi = -sphi^2/2 (== mean(phi) = 1,
 * matching the v.3 mphi convention).
 *
 * THREADING: the per-gene likelihood loop is partitioned via reduce_sum,
 * giving N-fold within-chain parallelism when sampling is invoked with
 * threads_per_chain > 1 (and the model is compiled with STAN_THREADS=true).
 * grainsize=1 (default) lets TBB choose chunks adaptively; explicit
 * grainsize ~ G / (2 * threads_per_chain) is sometimes faster for the
 * inner loop's relatively small per-gene cost.
 *
 * Priors (vector ops on dM, dEta, log_phi) stay outside reduce_sum -- Stan
 * already vectorizes those efficiently, and the per-gene non-iid work is
 * concentrated in the multinomial-logit normalization.
 *
 * Layout: per AA, the non-ref codons live in a contiguous slice of dM and
 * dEta from aa_start[a] to aa_end[a] (1-indexed inclusive).  Codon counts
 * y_k[g, k] are correspondingly stored in a flat per-gene * per-non-ref-codon
 * matrix.  Total AA-a residues per gene N_ga[g, a] are stored separately so
 * the normalization log_sum_exp(eta_full) can be subtracted exactly once per
 * (gene, AA) cell rather than once per codon.
 * ============================================================================ */

functions {
    // Per-gene multinomial-logit log-likelihood, summed over genes in the
    // slice that reduce_sum hands to this worker.  `slice_g` is an array of
    // 1-based gene IDs (a sub-range of the gene_indices array in transformed
    // data).  start/end are the original-array bounds and are unused here
    // because we iterate over slice_g directly.
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
            // Term 1 (vectorized): sum over non-ref codons of y_{g,k} * eta_k
            lp += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);
            // Term 2: per (g, a) cell normalization
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

    real<lower=0> sphi;

    int<lower=1> grainsize;                       // reduce_sum partition size
}

transformed data {
    // Gene-index array used as the slice variable for reduce_sum.
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] log_phi;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // Priors (vectorized, no benefit from threading)
    dM      ~ normal(dM_prior_mean, dM_prior_sd);
    dEta    ~ normal(dEta_prior_mean, dEta_prior_sd);
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);

    // Per-gene likelihood via reduce_sum (threaded if STAN_THREADS=true)
    target += reduce_sum(partial_sum_lpdf, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
