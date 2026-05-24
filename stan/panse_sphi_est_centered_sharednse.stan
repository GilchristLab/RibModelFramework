/* ============================================================================
 * panse_sphi_est_centered_sharednse.stan -- PANSE v2: log_phi sampled,
 * sphi ALSO sampled, single shared NSE scalar.  Centered parameterization.
 *
 * Promotes sphi from data (fixed) to a parameter relative to
 * panse_basic_sharednse.stan.  Diagnostic motivation: v1 (basic-sharednse)
 * showed slow joint mixing of phi[g] and lambdaPrime[c] due to the
 * multiplicative degeneracy mu = alpha * phi * sigma / (U * lambda).
 * Adding a hyperparameter on the phi distribution (sphi) lets the model
 * shift the multiplicative scale freely, which may improve the joint
 * geometry.  If v2 mixes cleanly, that tells us where in the model the
 * mixing problem actually lives (it's in the multiplicative scale, not
 * the per-codon parameters).
 *
 * Prior:
 *     sphi    ~ half-normal(0, sphi_prior_sd)     (real<lower=0>)
 *     log_phi ~ Normal(-sphi^2/2, sphi)          (v.3 mPhi convention,
 *                                                  mean(phi) = 1)
 *
 * See panse_basic_sharednse.stan for everything else; this file diffs by
 * 4 lines: (a) data declares sphi_prior_sd instead of sphi, (b) parameter
 * block adds real<lower=0> sphi, (c) model block adds the sphi prior.
 *
 * Noncentered companion (TODO): panse_sphi_est_noncentered_sharednse.stan,
 * for cases where the centered parameterization shows funnel-driven
 * divergences.
 * ============================================================================ */

functions {
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
    int<lower=1> G;
    int<lower=1> C;
    int<lower=G> P;

    array[G + 1] int<lower=1> gene_offset;
    array[P]     int<lower=1, upper=C> codon_at_pos;
    array[P]     int<lower=0> y;
    array[P]     int<lower=0, upper=1> like_mask;
    int<lower=0, upper=1> all_unmasked;

    real<lower=0> U;
    real<lower=0> sphi_prior_sd;                   // half-normal(0, .)

    real log_alpha_prior_mean;
    real<lower=0> log_alpha_prior_sd;
    real log_lambda_prior_mean;
    real<lower=0> log_lambda_prior_sd;

    real log_alpha_lower;
    real log_alpha_upper;
    real log_lambda_lower;
    real log_lambda_upper;

    real log_nse_lower;
    real log_nse_upper;
    int<lower=0, upper=1> nse_log_uniform;

    int<lower=0, upper=1> emit_log_lik;

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
    real log_U = log(U);
    int n_log_lik = emit_log_lik * P;
}

parameters {
    vector<lower=log_alpha_lower,  upper=log_alpha_upper>[C]  log_alpha;
    vector<lower=log_lambda_lower, upper=log_lambda_upper>[C] log_lambdaPrime;
    real<lower=log_nse_lower, upper=log_nse_upper> log_NSERate_shared;
    real<lower=0> sphi;
    vector[G] log_phi;
}

transformed parameters {
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    real NSERate_shared            = exp(log_NSERate_shared);
    vector<lower=0>[C] NSERate     = rep_vector(NSERate_shared, C);
    vector<lower=0>[G] phi         = exp(log_phi);

    vector[C] log_alpha_term;
    vector[C] log_psuccess;
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
    log_alpha       ~ normal(log_alpha_prior_mean,  log_alpha_prior_sd);
    log_lambdaPrime ~ normal(log_lambda_prior_mean, log_lambda_prior_sd);

    if (nse_log_uniform == 0) {
        target += log_NSERate_shared;
    }

    sphi    ~ normal(0, sphi_prior_sd);           // half-normal (lower=0)
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         gene_offset, codon_at_pos, y, like_mask, all_unmasked,
                         alpha, log_alpha_term, log_psuccess, log_phi);
}

generated quantities {
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
