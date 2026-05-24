/* ============================================================================
 * panse_csp_only.stan -- PANSE codon-specific parameters only (v0).
 *
 * Per-position negative-binomial-2 likelihood for ribosome footprint counts
 * along the ORF, with phi (gene synthesis rate) HELD FIXED as data.  Only
 * the codon-specific parameters Alpha[c], LambdaPrime[c], NSERate[c] are
 * sampled.  This is the validation gate before phi is added as a free
 * parameter in panse_basic.stan / panse_sphi_est_*.stan.
 *
 * Likelihood (per gene g, per position p with codon c = codon_at_pos[p]):
 *
 *     mu[g,p]  = alpha[c] * phi[g] * survive[g,p] / (U * lambdaPrime[c])
 *     y[g,p]  ~  NegBinomial2(mu[g,p], r = alpha[c])
 *
 * where survive[g,p] = P(ribosome reaches position p of gene g) is the
 * cumulative product of per-codon success probabilities at all upstream
 * positions p' < p.  In log space:
 *
 *     log_survive[g,1]   = 0
 *     log_survive[g,p+1] = log_survive[g,p] + log_psuccess(codon_at_pos[p])
 *     log_psuccess(c)    = -a_over_lv + a_over_lv / (lambdaPrime[c] * v)
 *                          + 0.5 * a_over_lv^2
 *     a_over_lv          = alpha[c] / (lambdaPrime[c] * v)
 *     v                  = 1 / NSERate[c]
 *
 * The log_psuccess form is a 2nd-order Taylor approximation that lives in
 * PANSEModel::elongationUntilIndexApproximation2ProbabilityLog (RMF C++).
 * In the C++, log_psuccess is clamped to 0 if positive (RMF safety net);
 * Stan omits that non-smooth clamp -- the posterior is expected to stay in
 * a regime where the approximation is well below 0.
 *
 * Partition function U:
 *     U = Z(mixture) / Y  is treated as fixed data here (computed once at
 *     data-prep time using the same formula as
 *     PANSEParameter::initializePartitionFunction).  Joint sampling of U
 *     is deferred to v3+.
 *
 * Data layout (CSR, NOT per-AA aggregation):
 *     positions live in one flat array y[P] of length P (sum of ORF lengths
 *     across all genes).  gene_offset[g] : gene_offset[g+1] - 1 gives the
 *     slice for gene g.  codon_at_pos[p] gives the 1-indexed codon ID at
 *     each position.  like_mask[p] = 1 means the position contributes to
 *     the likelihood; like_mask[p] = 0 means the position advances survive
 *     but does NOT contribute (RMF's "sigma-only" flag, positionMixture+1<0).
 *
 * NSE prior:
 *     log_NSERate is declared with bounds [log_nse_lower, log_nse_upper],
 *     giving an implicit Uniform prior on log_NSERate (== Log-Uniform on
 *     NSERate).  If nse_log_uniform == 0, an explicit Jacobian
 *     `target += sum(log_NSERate)` converts that to Natural-Uniform on
 *     NSERate.  The YAML key `nserate.prior.type` picks one of:
 *         Log-Uniform     -> nse_log_uniform = 1   (default)
 *         Natural-Uniform -> nse_log_uniform = 0
 *
 * PERFORMANCE: all per-codon work is precomputed once per draw in
 * transformed parameters (vectors of length C=61) rather than recomputed
 * inside the per-position inner loop (P ~ 747k positions for a full
 * genome).  This is ~5-10x faster than the naive version because it
 * eliminates ~12000x redundant log() / autodiff operations per iteration
 * and shrinks the autodiff tape correspondingly.
 *
 * THREADING: per-gene likelihood via reduce_sum.  Compile with
 *     cmdstan_model(..., cpp_options = list(stan_threads = TRUE))
 * and sample with threads_per_chain > 1.
 * ============================================================================ */

functions {
    /* Per-slice partial sum, summed over genes assigned to this worker.
     * slice_g is a sub-array of gene_indices.  Returns the partial
     * log-likelihood contribution.
     *
     * Per-codon precomputations (log_alpha_term, log_psuccess) are
     * hoisted to transformed parameters; per-gene work uses Stan's
     * vectorized neg_binomial_2_log_lpmf which fuses ~n NB2 evaluations
     * and their autodiff into a single call.  Together this is ~10-20x
     * faster than scalar-loop NB2 for genes with many positions.
     *
     * all_unmasked = 1: data has no sigma-only positions -- skip the
     *   per-position mask check; vectorize the whole gene at once.
     *   This is the common case (Weinberg, Wu, Mohammad CSVs all have
     *   like_mask == 1 for every position).
     * all_unmasked = 0: fall back to a scalar loop with the mask check.
     *   Only paid when actually needed. */
    real partial_sum(array[] int slice_g, int start, int end,
                     array[] int gene_offset,
                     array[] int codon_at_pos,
                     array[] int y,
                     array[] int like_mask,
                     int all_unmasked,
                     vector alpha,
                     vector log_alpha_term,
                     vector log_psuccess,
                     vector log_phi) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g  = slice_g[i];
            int p0 = gene_offset[g];
            int p1 = gene_offset[g + 1] - 1;
            int n  = p1 - p0 + 1;
            real lpg = log_phi[g];

            if (all_unmasked == 1) {
                /* Build per-gene log_mu and alpha vectors, then call
                 * vectorized NB2 once.  cumulative_sum on a shifted
                 * log_psuccess[] gives the per-position log_survive
                 * (with 0 prepended for position 1). */
                vector[n] log_psuccess_g;
                vector[n] alpha_g;
                for (j in 1:n) {
                    int c = codon_at_pos[p0 + j - 1];
                    log_psuccess_g[j] = log_psuccess[c];
                    alpha_g[j]        = alpha[c];
                }
                vector[n] log_survive_g;
                log_survive_g[1] = 0;
                if (n > 1)
                    log_survive_g[2:n] = cumulative_sum(log_psuccess_g[1:(n - 1)]);

                vector[n] log_mu_g;
                for (j in 1:n) {
                    int c = codon_at_pos[p0 + j - 1];
                    log_mu_g[j] = log_alpha_term[c] + lpg + log_survive_g[j];
                }
                lp += neg_binomial_2_log_lpmf(y[p0:p1] | log_mu_g, alpha_g);
            } else {
                /* Scalar loop with mask -- maintains log_survive across
                 * sigma-only positions.  Slower per-position but only
                 * used when the data actually contains masked positions. */
                real log_survive = 0;
                for (p in p0:p1) {
                    int c = codon_at_pos[p];
                    if (like_mask[p] == 1) {
                        lp += neg_binomial_2_log_lpmf(
                            y[p] | log_alpha_term[c] + lpg + log_survive,
                            alpha[c]);
                    }
                    log_survive += log_psuccess[c];
                }
            }
        }
        return lp;
    }
}

data {
    int<lower=1> G;                                // number of genes
    int<lower=1> C;                                // number of codons (61, no stops)
    int<lower=G> P;                                // total positions across all genes

    array[G + 1] int<lower=1> gene_offset;         // CSR offsets; gene g spans
                                                   //   gene_offset[g] : gene_offset[g+1]-1
    array[P]     int<lower=1, upper=C> codon_at_pos;
    array[P]     int<lower=0> y;                   // RFP counts
    array[P]     int<lower=0, upper=1> like_mask;  // 1 = include in likelihood
    int<lower=0, upper=1> all_unmasked;            // 1 = no sigma-only positions; enables
                                                   //   the per-gene vectorized NB2 fast path

    vector<lower=0>[G] phi;                        // synthesis rate per gene (FIXED, data)

    real<lower=0> U;                               // partition function (fixed)

    /* Priors on log-scale codon parameters.  YAML supplies natural-scale
     * lower/upper from RMF Uniform(0, 100); the R driver maps those to
     * (log_alpha_prior_mean, log_alpha_prior_sd) such that ~99% of the
     * weak-normal mass lands inside the natural-scale Uniform support. */
    real log_alpha_prior_mean;
    real<lower=0> log_alpha_prior_sd;
    real log_lambda_prior_mean;
    real<lower=0> log_lambda_prior_sd;

    /* Hard bounds on log_alpha / log_lambdaPrime to prevent HMC from
     * exploring extreme tails where neg_binomial_2_log_lpmf evaluates to
     * -nan during warmup proposals.  Set to log() of the YAML's natural-
     * scale uniform bounds; e.g. log(1e-3) ~ -6.9 to log(100) ~ 4.6. */
    real log_alpha_lower;
    real log_alpha_upper;
    real log_lambda_lower;
    real log_lambda_upper;

    /* NSE prior bounds (always declared in log space).  Matches RMF's
     * Log-Uniform(nserate.uniform.lower, nserate.uniform.upper). */
    real log_nse_lower;
    real log_nse_upper;
    int<lower=0, upper=1> nse_log_uniform;         // 1 = log-uniform, 0 = natural-uniform

    int<lower=0, upper=1> emit_log_lik;            // 1 = emit log_lik[P] in generated
                                                   //   quantities (24 GB at full Weinberg
                                                   //   scale); 0 = skip (codon-recovery
                                                   //   fits don't need WAIC).

    int<lower=1> grainsize;                        // reduce_sum partition size
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
    vector[G] log_phi = log(phi);
    real log_U = log(U);

    // Length of log_lik in generated quantities: P if emit_log_lik, else 0.
    int n_log_lik = emit_log_lik * P;
}

parameters {
    /* Log-transformed codon params with hard bounds at the natural-scale
     * uniform support (YAML alpha.prior.lower/upper).  Bounds prevent HMC
     * from exploring extreme tails where the NB2 likelihood evaluates to
     * -nan during warmup proposals; the weak normal prior still does the
     * within-bounds shaping. */
    vector<lower=log_alpha_lower,  upper=log_alpha_upper>[C]  log_alpha;
    vector<lower=log_lambda_lower, upper=log_lambda_upper>[C] log_lambdaPrime;

    /* NSE: declared in log space with bounds.  Implicit uniform prior in
     * log space gives Log-Uniform on NSERate; Jacobian below converts to
     * Natural-Uniform when nse_log_uniform == 0. */
    vector<lower=log_nse_lower, upper=log_nse_upper>[C] log_NSERate;
}

transformed parameters {
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    vector<lower=0>[C] NSERate     = exp(log_NSERate);

    /* Per-codon precomputations: hoisted out of the per-position inner
     * loop in partial_sum to avoid redundant log() / arithmetic calls
     * across all P positions (P ~ 747k for full Weinberg fit). */
    vector[C] log_alpha_term;       // = log_alpha[c] - log_U - log_lambdaPrime[c]
    vector[C] log_psuccess;         // 2nd-order Taylor of log P(success at codon c)
    for (c in 1:C) {
        log_alpha_term[c] = log_alpha[c] - log_U - log_lambdaPrime[c];
        real v         = 1.0 / NSERate[c];
        real a_over_lv = alpha[c] / (lambdaPrime[c] * v);
        log_psuccess[c] = -a_over_lv
                          + a_over_lv / (lambdaPrime[c] * v)
                          + 0.5 * a_over_lv * a_over_lv;
    }
}

model {
    /* Codon-level priors (vectorized, no benefit from threading) */
    log_alpha       ~ normal(log_alpha_prior_mean,  log_alpha_prior_sd);
    log_lambdaPrime ~ normal(log_lambda_prior_mean, log_lambda_prior_sd);

    /* NSE prior: implicit uniform on log_NSERate (= Log-Uniform on
     * NSERate).  Convert to Natural-Uniform when nse_log_uniform == 0. */
    if (nse_log_uniform == 0) {
        target += sum(log_NSERate);  // Jacobian: log|d NSERate / d log_NSERate|
    }

    /* Per-gene likelihood via reduce_sum (threaded if STAN_THREADS=true) */
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         gene_offset, codon_at_pos, y, like_mask, all_unmasked,
                         alpha, log_alpha_term, log_psuccess, log_phi);
}

generated quantities {
    /* Per-position log-likelihood for WAIC / LOO post-hoc.  Positions with
     * like_mask == 0 contribute 0 (no penalty).
     * Emitted only when emit_log_lik == 1; at full Weinberg scale
     * (P ~ 747k, 4000 draws) the GQ output is ~24 GB and the per-sample
     * GQ pass adds ~30% wall time. */
    vector[n_log_lik] log_lik;
    if (emit_log_lik) {
        for (g in 1:G) {
            int p0 = gene_offset[g];
            int p1 = gene_offset[g + 1] - 1;
            real log_survive = 0;
            real lpg = log_phi[g];
            for (p in p0:p1) {
                int c = codon_at_pos[p];
                if (like_mask[p] == 1) {
                    log_lik[p] = neg_binomial_2_log_lpmf(
                        y[p] | log_alpha_term[c] + lpg + log_survive,
                        alpha[c]);
                } else {
                    log_lik[p] = 0;
                }
                log_survive += log_psuccess[c];
            }
        }
    }
}
