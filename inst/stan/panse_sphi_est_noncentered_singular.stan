/* ============================================================================
 * panse_sphi_est_noncentered_singular.stan -- SINGULAR (deterministic) dwell
 * variant of panse_sphi_est_noncentered.stan.
 *
 * DWELL DISTRIBUTION: the per-codon dwell time W is taken DETERMINISTIC at its
 * mean W = E[W] = alpha/lambda.  This is the alpha->inf (zero dwell-variance)
 * limit of the Gamma-dwell model, with the Gamma mean held fixed so the
 * parameters keep their meaning.  Two consequences vs the Gamma parent:
 *
 *   (1) SURVIVAL.  With W deterministic, E_W[v/(W+v)] is the integrand at the
 *       point value:  psuccess = 1/(1+q),  q = alpha*NSERate/lambda,
 *       log_psuccess = -log1p(q).  (Same as the pointval surrogate.)
 *
 *   (2) COUNTS.  With W deterministic there is nothing to marginalize:
 *       Y ~ Poisson(mu),  mu = phi*alpha/(U*lambda)*survival  (equidispersed).
 *       This REPLACES the parent's NegBin2(mu, size=alpha).
 *
 * So this model = Poisson counts + POINT-VALUE survival, the self-consistent
 * deterministic-dwell model.  (It is NOT panse_nc_poisson.stan, which mixes
 * Poisson counts with the random-dwell Gamma survival.)
 *
 * The parameter vector is byte-identical to panse_sphi_est_noncentered.stan
 * (log_alpha[C], log_lambdaPrime[C], log_NSERate[C], sphi, z_phi[G]); the
 * size=alpha argument of the count lpmf is simply dropped (Poisson has no
 * dispersion parameter).  Everything else (noncentered log_phi, generalized
 * phi prior, soft mean(phi)=1 anchor, partial_sum, generated quantities) is
 * identical to the parent.  See notes/dwell-time-distributions.md.
 * ============================================================================ */

functions {
    // Point-value (deterministic-dwell) log survival: E_W[v/(W+v)] at W=E[W].
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
                for (j in 1:n) {
                    int c = codon_at_pos[p0 + j - 1];
                    log_psuccess_g[j] = log_psuccess[c];
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
                lp += poisson_log_lpmf(y[p0:p1] | log_mu_g);   // deterministic dwell -> Poisson
            } else {
                real log_survive = 0;
                for (p in p0:p1) {
                    int c = codon_at_pos[p];
                    if (like_mask[p] == 1) {
                        lp += poisson_log_lpmf(
                            y[p] | log_alpha_term[c] + lpg + log_survive);
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
        // Deterministic-dwell (point-value) survival: 1/(1+q).
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
                    log_lik[p] = poisson_log_lpmf(
                        y[p] | log_alpha_term[c] + lpg + log_survive);
                } else {
                    log_lik[p] = 0;
                }
                log_survive += log_psuccess[c];
            }
        }
    }
}
