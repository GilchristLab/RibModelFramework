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
 * THREADING: per-gene likelihood via reduce_sum.  Compile with
 *     cmdstan_model(..., cpp_options = list(stan_threads = TRUE))
 * and sample with threads_per_chain > 1.
 * ============================================================================ */

functions {
    /* Per-slice partial sum, summed over genes assigned to this worker.
     * slice_g is a sub-array of gene_indices.  Returns the partial
     * log-likelihood contribution. */
    real partial_sum(array[] int slice_g, int start, int end,
                     array[] int gene_offset,
                     array[] int codon_at_pos,
                     array[] int y,
                     array[] int like_mask,
                     vector alpha, vector lambdaPrime, vector NSERate,
                     vector phi, vector log_phi,
                     real U, real log_U) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g  = slice_g[i];
            int p0 = gene_offset[g];
            int p1 = gene_offset[g + 1] - 1;
            real log_survive = 0;  // log P(reach position p0); position 1 has prob 1.

            for (p in p0:p1) {
                int  c   = codon_at_pos[p];
                real lpd = lambdaPrime[c];
                real a   = alpha[c];
                real nse = NSERate[c];

                if (like_mask[p] == 1) {
                    /* mu = a * phi[g] * exp(log_survive) / (U * lpd)
                     * Stan's neg_binomial_2_log_lpmf takes log_mu and is more
                     * numerically stable when log_survive is very negative. */
                    real log_mu = log(a) + log_phi[g] + log_survive
                                  - log_U - log(lpd);
                    lp += neg_binomial_2_log_lpmf(y[p] | log_mu, a);
                }

                /* Update log_survive for the NEXT position (after the
                 * likelihood evaluation -- matches RMF's off-by-one). */
                real v         = 1.0 / nse;
                real a_over_lv = a / (lpd * v);
                log_survive   += -a_over_lv
                                 + a_over_lv / (lpd * v)
                                 + 0.5 * a_over_lv * a_over_lv;
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

    /* NSE prior bounds (always declared in log space).  Matches RMF's
     * Log-Uniform(nserate.uniform.lower, nserate.uniform.upper). */
    real log_nse_lower;
    real log_nse_upper;
    int<lower=0, upper=1> nse_log_uniform;         // 1 = log-uniform, 0 = natural-uniform

    int<lower=1> grainsize;                        // reduce_sum partition size
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
    vector[G] log_phi = log(phi);
    real log_U = log(U);
}

parameters {
    /* Log-transformed codon params: weak normal prior in log space, no
     * upper bound (RMF's Uniform(0, 100) is matched loosely via the prior
     * SDs supplied by the YAML).  Log-transform gives HMC smooth
     * gradients across many orders of magnitude. */
    vector[C] log_alpha;
    vector[C] log_lambdaPrime;

    /* NSE: declared in log space with bounds.  Implicit uniform prior in
     * log space gives Log-Uniform on NSERate; Jacobian below converts to
     * Natural-Uniform when nse_log_uniform == 0. */
    vector<lower=log_nse_lower, upper=log_nse_upper>[C] log_NSERate;
}

transformed parameters {
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    vector<lower=0>[C] NSERate     = exp(log_NSERate);
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
                         gene_offset, codon_at_pos, y, like_mask,
                         alpha, lambdaPrime, NSERate,
                         phi, log_phi,
                         U, log_U);
}

generated quantities {
    /* Per-position log-likelihood for WAIC / LOO post-hoc.  Positions with
     * like_mask == 0 contribute 0 (no penalty). */
    vector[P] log_lik;
    {
        for (g in 1:G) {
            int p0 = gene_offset[g];
            int p1 = gene_offset[g + 1] - 1;
            real log_survive = 0;
            for (p in p0:p1) {
                int  c   = codon_at_pos[p];
                real lpd = lambdaPrime[c];
                real a   = alpha[c];
                real nse = NSERate[c];

                if (like_mask[p] == 1) {
                    real log_mu = log(a) + log_phi[g] + log_survive
                                  - log_U - log(lpd);
                    log_lik[p] = neg_binomial_2_log_lpmf(y[p] | log_mu, a);
                } else {
                    log_lik[p] = 0;
                }

                real v         = 1.0 / nse;
                real a_over_lv = a / (lpd * v);
                log_survive   += -a_over_lv
                                 + a_over_lv / (lpd * v)
                                 + 0.5 * a_over_lv * a_over_lv;
            }
        }
    }
}
