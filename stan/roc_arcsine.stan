/* ============================================================================
 * roc_arcsine.stan  --  ROC model with hybrid arcsine / exact-multinomial
 *                       likelihood.
 *
 * For each (gene g, amino acid group a) pair:
 *
 *   N = N_ga[g, a]  (total codon count -- fixed data, never a parameter)
 *
 *   N == 0               -> skip
 *   N >= approx_min_n    -> arcsine approximation (K_a - 1 marginal binomials)
 *   N <  approx_min_n    -> exact multinomial
 *
 * The threshold is on fixed data so no gradient discontinuity exists for HMC.
 * The default of approx_min_n = 20 matches the C++ implementation (src/ROCModel.cpp).
 *
 * Arcsine log-likelihood for one (gene, AA) pair:
 *
 *   logL = sum_{k=non-ref} -2N * (asin(sqrt(c_k/N)) - asin(sqrt(p_k)))^2
 *
 * Only K_a - 1 non-reference codons are summed.  The reference codon is
 * excluded because asin(sqrt(x)) + asin(sqrt(1-x)) = pi/2 identically,
 * so the K-th arcsine term is a linear function of the first K-1 terms --
 * including it would double-count the constraint sum(c_k) = N.
 *
 * Exact log-likelihood (used when N < approx_min_n):
 *
 *   logL = sum_k c_k * log(p_k)
 *        = sum_{k=non-ref} c_k * eta_k - N * log_sum_exp(eta_full)
 *
 * where eta_k = -dM[k] - dEta[k]*phi[g]  and  eta_ref = 0.
 *
 * Data layout (same as roc_sphi_est.stan):
 *   y_k[G, K]    -- non-reference codon counts; K = total non-ref codons
 *   N_ga[G, A]   -- total codon count per gene per AA group (all codons incl. ref)
 *   aa_start[A], aa_end[A]  -- inclusive index range into y_k / dM / dEta
 *
 * Parametrization options (same flags as roc_sphi_est.stan):
 *   noncentered = 0  centered    (default; best when data strongly anchors phi)
 *   noncentered = 1  non-centered (better for data-sparse fits; removes funnel)
 *   anchor_phi  = 0  mean(phi) = 1 via mphi = -sphi^2/2  (default)
 *   anchor_phi  = 1  soft median(phi) ~ 1 via prior on mphi_param
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
 * threads_per_chain > 1 to enable reduce_sum parallelism.
 * ============================================================================ */

functions {

    /* partial_sum_hybrid: per-gene log-likelihood contribution.
     *
     * For each gene in slice_g, iterates over AA groups and applies
     * either the arcsine approximation (N >= approx_min_n) or the
     * exact multinomial (N < approx_min_n).
     */
    real partial_sum_hybrid(
            array[] int slice_g, int start, int end,
            int A,
            array[] int aa_start,
            array[] int aa_end,
            array[,] int y_k,
            array[,] int N_ga,
            vector dM,
            vector dEta,
            vector phi,
            int approx_min_n,
            array[,] real asin_sqrt_phat) {

        real lp = 0;

        for (i in 1:size(slice_g)) {
            int g = slice_g[i];

            for (a in 1:A) {
                int N = N_ga[g, a];
                if (N == 0) continue;

                int s         = aa_start[a];
                int e         = aa_end[a];
                int n_nonref  = e - s + 1;   // K_a - 1 (non-reference codons)

                /* Build unnormalized log-numerator vector.
                 * Index 1 = reference codon (eta = 0).
                 * Indices 2..n_nonref+1 = non-reference codons (eta = -dM - dEta*phi). */
                vector[n_nonref + 1] eta_full;
                eta_full[1] = 0.0;
                for (k in 1:n_nonref)
                    eta_full[k + 1] = -dM[s - 1 + k] - dEta[s - 1 + k] * phi[g];

                real log_Z = log_sum_exp(eta_full);

                if (N >= approx_min_n) {
                    /* ---- Arcsine approximation --------------------------------
                     * Sum K_a - 1 marginal binomial arcsine terms.
                     * Reference codon excluded; see file header for rationale.
                     *
                     * p_k is clamped to (1e-12, 1-1e-12) to prevent a gradient
                     * singularity.  The cross-gradient d/d(eta_j)[arcsine_k term]
                     * for j != k is proportional to p_j * sqrt(p_k) / sqrt(1-p_k),
                     * which diverges as p_k -> 1.  The diagonal gradient cancels
                     * via the softmax Jacobian and is always finite, but the cross
                     * terms do not.  Clamping at 1e-12 keeps gradients bounded
                     * everywhere without affecting inference at reasonable parameter
                     * values (p_k = 1-1e-12 corresponds to a logit of ~27).
                     *
                     * asin_sqrt_phat[g, s-1+k] = asin(sqrt(c_k/N)) is precomputed
                     * in transformed data -- it depends only on fixed data. */
                    real N_r = N;
                    for (k in 1:n_nonref) {
                        real p_k  = fmin(1.0 - 1e-12, fmax(1e-12,
                                         exp(eta_full[k + 1] - log_Z)));
                        real diff = asin_sqrt_phat[g, s - 1 + k] - asin(sqrt(p_k));
                        lp += -2.0 * N_r * diff ^ 2;
                    }
                } else {
                    /* ---- Exact multinomial ------------------------------------
                     * logL = sum_{k=non-ref} c_k * eta_k  -  N * log_Z
                     * The -N*log_Z term includes the reference codon's -c_ref*log_Z
                     * contribution since N = sum of ALL K_a codons. */
                    for (k in 1:n_nonref)
                        lp += y_k[g, s - 1 + k] * eta_full[k + 1];
                    lp -= N * log_Z;
                }
            }
        }
        return lp;
    }
}

data {
    int<lower=1> G;                        // number of genes
    int<lower=1> A;                        // number of AA groups (synonymous families)
    int<lower=1> K;                        // total non-reference codons across all AAs
    array[A] int<lower=1> aa_start;        // first index into dM/dEta/y_k for each AA
    array[A] int<lower=1> aa_end;          // last  index into dM/dEta/y_k for each AA
    array[G, K] int<lower=0> y_k;          // non-reference codon counts [gene, codon]
    array[G, A] int<lower=0> N_ga;         // total codon count [gene, AA group]

    int<lower=1>     approx_min_n;         // min N for arcsine (default 20)

    vector[K]           dM_prior_mean;
    vector<lower=0>[K]  dM_prior_sd;
    vector[K]           dEta_prior_mean;
    vector<lower=0>[K]  dEta_prior_sd;
    real          sphi_prior_mean;
    real<lower=0> sphi_prior_sd;

    int<lower=0, upper=1> noncentered;     // 0 = centered (default), 1 = non-centered
    int<lower=0, upper=1> anchor_phi;      // 0 = mean(phi)=1, 1 = soft median(phi)~1
    real<lower=0> mphi_prior_sd;           // SD for mphi soft prior; ignored when anchor_phi=0
    int<lower=1>  grainsize;               // reduce_sum grain size
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;

    // Precompute asin(sqrt(c_k / N)) -- depends only on fixed data.
    // Eliminates one asin + sqrt per (gene, non-ref codon) per gradient call.
    array[G, K] real asin_sqrt_phat;
    for (g in 1:G) {
        for (a in 1:A) {
            int N = N_ga[g, a];
            int s = aa_start[a];
            int e = aa_end[a];
            for (k in 1:(e - s + 1)) {
                asin_sqrt_phat[g, s - 1 + k] =
                    (N == 0) ? 0.0 : asin(sqrt(y_k[g, s - 1 + k] * 1.0 / N));
            }
        }
    }
}

parameters {
    vector[K]   dM;
    vector[K]   dEta;
    vector[G]   latent_phi;               // log_phi (centered) or z_phi (non-centered)
    real<lower=0> sphi;
    real mphi_param;                      // used only when anchor_phi=1; phantom otherwise
}

transformed parameters {
    real   mphi    = (anchor_phi == 1) ? mphi_param : -0.5 * sphi * sphi;
    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    // ---- Priors -------------------------------------------------------
    dM   ~ normal(dM_prior_mean,   dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi ~ normal(sphi_prior_mean, sphi_prior_sd);   // positive via lower=0 constraint

    // mphi_param: active when anchor_phi=1; phantom (wide prior) otherwise
    if (anchor_phi == 1)
        mphi_param ~ normal(0, mphi_prior_sd);
    else
        mphi_param ~ normal(0, 1);

    // phi prior
    if (noncentered)
        latent_phi ~ std_normal();
    else
        latent_phi ~ normal(mphi, sphi);

    // ---- Hybrid arcsine / exact likelihood ----------------------------
    target += reduce_sum(partial_sum_hybrid, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi, approx_min_n, asin_sqrt_phat);
}

generated quantities {
    vector[G] phi_out = exp(log_phi);     // synthesis rates on original scale
}
