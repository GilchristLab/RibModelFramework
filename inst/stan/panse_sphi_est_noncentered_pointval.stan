/* ============================================================================
 * panse_sphi_est_noncentered_pointval.stan -- point-value-dwell SURROGATE of
 * panse_sphi_est_noncentered.stan.
 *
 * PURPOSE: a geometrically near-identical, cheap-to-fit surrogate whose
 * converged HMC posterior covariance seeds the dense_e inv_metric (and init)
 * of the full Gamma-dwell model.  The parameter vector is BYTE-IDENTICAL to
 * panse_sphi_est_noncentered.stan (same blocks, same ordering:
 * log_alpha[C], log_lambdaPrime[C], log_NSERate[C], sphi, z_phi[G]), so the
 * unconstrained posterior covariance transplants 1:1 -- no expand heuristic.
 *
 * ONLY difference from the Gamma model: the per-codon survival probability.
 * Instead of integrating v/(W+v) over W ~ Gamma(shape=alpha, rate=lambda)
 * (the upper-incomplete-gamma / Lentz continued fraction in the parent model),
 * the dwell is taken at its point value w = E[W] = alpha/lambda, moving the
 * expectation inside:
 *     psuccess_c = v / (w_c + v) = 1 / (1 + q_c),  q_c = alpha_c*NSERate_c/lambda_c
 *     log_psuccess_c = -log1p(q_c)
 * This is bounded in (0,1], monotone decreasing in NSERate, and -- unlike the
 * 2nd-order Taylor it ultimately descends from -- CANNOT breach psuccess>1.
 * It is the integrand evaluated at W=E[W], so it is the closest cheap collapse
 * of the Gamma survival (agrees to first order in q; q = alpha*NSERate/lambda
 * is << 1 in the relevant regime).
 *
 * The count model is UNCHANGED: NegBin2(mu, size=alpha) with the same
 * log_alpha_term mean.  We deliberately keep size=alpha so the surrogate shares
 * the Gamma model's parameter identifiability exactly; the surrogate is a
 * metric/init source, not a competing scientific model.
 *
 * Everything else (noncentered log_phi, generalized phi prior, soft mean(phi)=1
 * anchor, partial_sum, generated quantities) is identical to the parent.
 * ============================================================================ */

functions {
    // Point-value log survival: log E_W[v/(W+v)] approximated at W=E[W]=alpha/lambda.
    //   psuccess = 1/(1+q),  q = alpha*NSERate/lambda  ->  log_psuccess = -log1p(q)
    real log_psuccess_pointval(real alpha, real lambda, real NSERate) {
        return -log1p(alpha * NSERate / lambda);
    }

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
    real<lower=0> sphi_prior_sd;

    // Generalized phi prior.  phi_use_data=0: sphi estimated, per-population
    // hierarchical prior (original model).  phi_use_data=1: per-gene prior
    // mean and SD from data vectors; sphi still sampled but does not affect
    // log_phi computation.
    int<lower=0, upper=1> phi_use_data;
    vector[G] phi_prior_mu;      // per-gene prior mean  (used when phi_use_data=1)
    vector[G] phi_prior_sigma;   // per-gene prior SD    (used when phi_use_data=1)

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
    vector<lower=log_nse_lower, upper=log_nse_upper>[C] log_NSERate;
    real<lower=0> sphi;
    vector[G] z_phi;                                // noncentered: unit normals
}

transformed parameters {
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    vector<lower=0>[C] NSERate     = exp(log_NSERate);
    // Noncentered log_phi.
    // phi_use_data=0: population prior; log_phi ~ Normal(-0.5*sphi^2, sphi).
    // phi_use_data=1: per-gene prior;   log_phi[g] ~ Normal(phi_prior_mu[g], phi_prior_sigma[g]).
    vector[G] log_phi;
    if (phi_use_data) {
        log_phi = phi_prior_mu + phi_prior_sigma .* z_phi;
    } else {
        log_phi = -0.5 * sphi * sphi + sphi * z_phi;
    }
    vector<lower=0>[G] phi = exp(log_phi);

    vector[C] log_alpha_term;
    vector[C] log_psuccess;
    for (c in 1:C) {
        log_alpha_term[c] = log_alpha[c] - log_U - log_lambdaPrime[c];
        // Point-value survival surrogate (1/(1+q)); replaces the Gamma-integrated
        // upper-incomplete-gamma form of the parent model.  See header.
        log_psuccess[c] = log_psuccess_pointval(alpha[c], lambdaPrime[c], NSERate[c]);
    }
}

model {
    log_alpha       ~ normal(log_alpha_prior_mean,  log_alpha_prior_sd);
    log_lambdaPrime ~ normal(log_lambda_prior_mean, log_lambda_prior_sd);

    if (nse_log_uniform == 0) {
        target += sum(log_NSERate);
    }

    sphi  ~ normal(0, sphi_prior_sd);
    z_phi ~ std_normal();

    /* Soft mean(phi)=1 anchor; breaks the phi <-> lambda multiplicative ridge.
     * Applied only in the population (phi_use_data=0) parameterization, where
     * log_phi = -0.5*sphi^2 + sphi*z_phi and the anchor constrains mean(z_phi)~0.
     * In the data parameterization (phi_use_data=1), phi_prior_mu pins the mean. */
    if (!phi_use_data) {
        target += -0.5 * square((mean(log_phi) + 0.5 * sphi * sphi) / 0.01);
    }

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
