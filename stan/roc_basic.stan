/* ============================================================================
 * roc_basic.stan -- ROC model in Stan, first cut.
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
 * (not free parameters).  sphi is held fixed for this first cut, with the
 * lognormal phi prior constraint mphi = -sphi^2/2 (== mean(phi) = 1,
 * matching the v.3 mphi convention).
 *
 * Priors are simple Gaussians; per-codon mean and SD vectors are passed in
 * via data so we can match any v.5 YAML prior config (gcBias, etc.) without
 * re-implementing those choices inside the Stan model.
 *
 * Layout: per AA, the non-ref codons live in a contiguous slice of dM and
 * dEta from aa_start[a] to aa_end[a] (1-indexed inclusive).  Codon counts
 * y_k[g, k] are correspondingly stored in a flat per-gene * per-non-ref-codon
 * matrix.  Total AA-a residues per gene N_ga[g, a] are stored separately so
 * the normalization log_sum_exp(eta_full) can be subtracted exactly once per
 * (gene, AA) cell rather than once per codon.
 * ============================================================================ */

data {
    int<lower=1> G;                              // number of genes
    int<lower=1> A;                              // number of AAs (split.serine -> 22 typical)
    int<lower=1> K;                              // sum of (n_codons - 1) across all AAs

    // Ragged AA layout: for AA a, non-ref codons sit at indices aa_start[a]..aa_end[a]
    // (1-indexed, inclusive) within the dM/dEta vectors.
    array[A] int<lower=1> aa_start;
    array[A] int<lower=1> aa_end;

    // Codon counts.  y_k[g, k] is the count of non-ref codon k in gene g.
    // N_ga[g, a] is the total AA-a residue count in gene g (ref + non-ref).
    array[G, K] int<lower=0> y_k;
    array[G, A] int<lower=0> N_ga;

    // Priors (per-codon mean + SD, per the v.5 YAML schema).
    vector[K] dM_prior_mean;
    vector<lower=0>[K] dM_prior_sd;
    vector[K] dEta_prior_mean;
    vector<lower=0>[K] dEta_prior_sd;

    real<lower=0> sphi;                          // fixed; will be estimated in later phases
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
    // Priors
    dM ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Phi prior: lognormal with mean(phi)=1 -> mean(log phi) = -sphi^2/2
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);

    // Likelihood: contribute log P(y_{g,a,:}) for each (gene, AA) cell.
    // Decomposed as: sum_k y_k * eta_k  -  sum_{(g,a)} N_ga * log_sum_exp(eta_full)
    for (g in 1:G) {
        // Term 1 (vectorized): sum over non-ref codons of y_{g,k} * (-dM_k - dEta_k * phi_g)
        target += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);

        // Term 2: per (g, a) cell normalization.  Skip cells with no observations.
        for (a in 1:A) {
            if (N_ga[g, a] == 0) continue;
            int s = aa_start[a];
            int e = aa_end[a];
            int n = e - s + 1;
            vector[n + 1] eta_full;
            eta_full[1] = 0;                     // reference codon
            for (k in 1:n) {
                eta_full[k + 1] = -dM[s - 1 + k] - dEta[s - 1 + k] * phi[g];
            }
            target += -N_ga[g, a] * log_sum_exp(eta_full);
        }
    }
}
