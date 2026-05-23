/* ============================================================================
 * panse_basic.stan -- PANSE v1: log_phi sampled, sphi fixed.
 *
 * Extends panse_csp_only.stan by promoting phi from data to a parameter.
 * Prior:
 *     log_phi ~ Normal(-0.5 * sphi^2, sphi)    (lognormal mean(phi) = 1,
 *                                               v.3 mPhi convention)
 * sphi is held fixed as data.  See panse_sphi_est_centered.stan for the
 * variant that estimates sphi.
 *
 * See panse_csp_only.stan for full math derivation, data layout, and
 * NSE prior switch.
 * ============================================================================ */

functions {
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
            real log_survive = 0;

            for (p in p0:p1) {
                int  c   = codon_at_pos[p];
                real lpd = lambdaPrime[c];
                real a   = alpha[c];
                real nse = NSERate[c];

                if (like_mask[p] == 1) {
                    real log_mu = log(a) + log_phi[g] + log_survive
                                  - log_U - log(lpd);
                    lp += neg_binomial_2_log_lpmf(y[p] | log_mu, a);
                }

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
    int<lower=1> G;
    int<lower=1> C;
    int<lower=G> P;

    array[G + 1] int<lower=1> gene_offset;
    array[P]     int<lower=1, upper=C> codon_at_pos;
    array[P]     int<lower=0> y;
    array[P]     int<lower=0, upper=1> like_mask;

    /* phi removed from data; log_phi is now a parameter (see below). */

    real<lower=0> U;
    real<lower=0> sphi;                            // synthesis-rate stddev (fixed)

    real log_alpha_prior_mean;
    real<lower=0> log_alpha_prior_sd;
    real log_lambda_prior_mean;
    real<lower=0> log_lambda_prior_sd;

    real log_nse_lower;
    real log_nse_upper;
    int<lower=0, upper=1> nse_log_uniform;

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
    real log_U = log(U);
}

parameters {
    vector[C] log_alpha;
    vector[C] log_lambdaPrime;
    vector<lower=log_nse_lower, upper=log_nse_upper>[C] log_NSERate;
    vector[G] log_phi;
}

transformed parameters {
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    vector<lower=0>[C] NSERate     = exp(log_NSERate);
    vector<lower=0>[G] phi         = exp(log_phi);
}

model {
    log_alpha       ~ normal(log_alpha_prior_mean,  log_alpha_prior_sd);
    log_lambdaPrime ~ normal(log_lambda_prior_mean, log_lambda_prior_sd);

    if (nse_log_uniform == 0) {
        target += sum(log_NSERate);
    }

    /* lognormal phi prior with mean = 1 (v.3 mPhi convention) */
    log_phi ~ normal(-0.5 * sphi * sphi, sphi);

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         gene_offset, codon_at_pos, y, like_mask,
                         alpha, lambdaPrime, NSERate,
                         phi, log_phi,
                         U, log_U);
}

generated quantities {
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
