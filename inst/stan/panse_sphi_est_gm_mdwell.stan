/* ============================================================================
 * panse_sphi_est_gm_mdwell.stan -- two stacked reparameterizations aimed at the
 *   PANSE posterior ridges: a GEOMETRIC-mean-one sum-to-zero phi scale-pin
 *   (inherited from panse_sphi_est_sumzero_geomean.stan) PLUS a MEAN-DWELL
 *   rotation of the per-codon (alpha, lambda) basis.
 *
 * MEAN-DWELL ROTATION: instead of sampling log_lambdaPrime[C] directly, sample
 *   log_mdwell = log(alpha/lambda) (the per-codon mean dwell time) and derive
 *   log_lambdaPrime = log_alpha - log_mdwell.  alpha and lambda are strongly
 *   correlated along the mean-dwell direction (the alpha<->lambda ridge, a
 *   documented PANSE mixing bottleneck); sampling in (log_alpha, log_mdwell)
 *   turns that direction into a coordinate axis, decorrelating it.  This
 *   rotation is not present in any other model in the family.
 *
 * STATUS: EXPERIMENTAL, UNTESTED.  Not registered in .PANSE_STAN_MODELS, so
 *   fit_panse_stan() cannot dispatch it yet -- wire + smoke-test before any
 *   use.  Generated 2026-06-07 during the scale-gauge/coupling exploration;
 *   committed for preservation, not validated.
 *
 * --- geomean scale-pin (inherited) -------------------------------------------
 * MOTIVATION (2026-06-06 coupling analysis): the parent sumzero model pins
 * mean(log_phi) = -0.5*sphi^2, so sphi shifts the phi LOCATION; the global
 * lambda scale must track that shift (Z ~ sum alpha/lambda), producing a
 * dominant sphi<->lambda posterior ridge (mean|cor| ~ 0.80) that keeps NUTS
 * treedepth-saturated.  This variant drops the -0.5*sphi^2 offset so
 * mean(log_phi) = 0 regardless of sphi (sphi = pure spread), aiming to break
 * that ridge.  Absolute phi scale is a gauge (likelihood sees phi/U), so
 * arithmetic-mean-one -> geometric-mean-one is scientifically benign; only the
 * normalization convention changes (E[phi]=exp(0.5 sphi^2) here vs 1 in parent).
 *
 * Original sumzero header follows.
 * ----------------------------------------------------------------------------
 * panse_sphi_est_sumzero.stan -- per-codon NSE, SUM-TO-ZERO phi scale-pin.
 *
 * HARD scale-pin variant of panse_sphi_est_noncentered.stan, and the Stan
 * analog of the native MCMC's fix.z.  The PANSE likelihood depends only on
 * phi_g / U, so the overall phi scale is a likelihood-flat gauge (the partition
 * function Z is a pure normalization).  The noncentered model leaves this gauge
 * free -- mean(z_phi) drifts and couples to sphi (the sphi<->mean_log_phi ridge,
 * cor ~ -0.98, the dominant mixing bottleneck) -- and tries to tame it with a
 * SOFT anchor (mean(log_phi) penalty, SD 0.01), which is a footgun: too loose
 * and the scale still floats, too tight and it injects high curvature HMC chokes
 * on.
 *
 * This model removes the gauge BY CONSTRUCTION instead: z_phi is a
 * sum_to_zero_vector, so mean(z_phi) = 0 exactly, mean(log_phi) is pinned to
 * -0.5*sphi^2 (the prior-consistent scale, geometric-mean phi = exp(-0.5 sphi^2)),
 * and the flat direction is deleted from the geometry -- no penalty curvature.
 * This is the reparameterization (hard) equivalent of fixing Z; the soft anchor
 * is dropped entirely.  Known minor cost: the sum_to_zero (G-1)/G variance
 * reduction gives a small sphi upward bias at very low G (~0.56 vs 0.50 at
 * G=100; negligible at production G).
 *
 * Differs from panse_sphi_est_noncentered.stan ONLY by:
 *   parameters:  vector[G] z_phi  ->  sum_to_zero_vector[G] z_phi
 *   model:       drop the soft mean(phi)=1 anchor term
 * Per-codon NSE, Lentz-CF survival, generalized phi prior all identical.
 * ============================================================================ */

functions {
    // log of the UNregularized upper incomplete gamma Gamma(s, x), x > 0, via
    // the forward modified-Lentz continued fraction (see noncentered parent).
    real log_upper_incomplete_gamma(real s, real x) {
        real tol   = 1e-10;
        real FPMIN = 1e-300;
        real b = x + 1 - s;
        if (abs(b) < FPMIN) b = FPMIN;
        real c = 1 / FPMIN;
        real d = 1 / b;
        real h = d;
        int  i = 1;
        real del = 0;
        int  converged = 0;
        while (i <= 2000 && converged == 0) {
            real an = -i * (i - s);
            b = b + 2;
            d = an * d + b;  if (abs(d) < FPMIN) d = FPMIN;
            c = b + an / c;  if (abs(c) < FPMIN) c = FPMIN;
            d = 1 / d;
            del = d * c;
            h = h * del;
            if (abs(del - 1) < tol) converged = 1;
            i = i + 1;
        }
        if (converged == 0)
            reject("log_upper_incomplete_gamma: CF did not converge for s=", s, " x=", x);
        return s * log(x) - x + log(h);
    }

    real log_psuccess_hybrid(real alpha, real lambda, real NSERate) {
        real v  = 1.0 / NSERate;
        real lv = lambda * v;
        real q  = alpha * NSERate / lambda;
        if (q < 0.005) {
            real a_over_lv = alpha / lv;
            return -a_over_lv + a_over_lv / lv + 0.5 * a_over_lv * a_over_lv;
        }
        return alpha * log(lv) + lv + log_upper_incomplete_gamma(1 - alpha, lv);
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

    // Generalized phi prior.  phi_use_data=0: population hierarchical prior;
    // phi_use_data=1: per-gene prior mean/SD from data (sphi still sampled).
    int<lower=0, upper=1> phi_use_data;
    vector[G] phi_prior_mu;
    vector[G] phi_prior_sigma;

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
    // MEAN-DWELL rotation: sample log_mdwell = log(alpha/lambda) = log_alpha -
    // log_lambdaPrime instead of log_lambdaPrime.  The within-codon
    // (log_alpha,log_lambda) ridge is LINEAR with major axis ~48.8deg (the (1,1)
    // direction) and condition number ~17: alpha/lambda (= E[dwell]) is well
    // pinned, the common magnitude is the loose ridge.  In (log_alpha, log_mdwell)
    // coords the ridge is axis-aligned, so a DIAGONAL metric (diag_e) suffices --
    // no O(n_par^2) dense_e + ever-growing warmup as the genome grows.  The
    // log_lambda box bound is dropped (non-binding: lambda ~0.08-0.12, log ~-2.5,
    // well inside [-6.9,4.6]); the normal prior on the DERIVED log_lambdaPrime is
    // preserved, so the posterior is unchanged (pure reparameterization).
    vector[C] log_mdwell;
    vector<lower=log_nse_lower, upper=log_nse_upper>[C] log_NSERate;
    real<lower=0> sphi;
    sum_to_zero_vector[G] z_phi;   // mean(z_phi)=0 by construction -> scale pinned
}

transformed parameters {
    vector[C]          log_lambdaPrime = log_alpha - log_mdwell;   // ridge rotation
    vector<lower=0>[C] alpha       = exp(log_alpha);
    vector<lower=0>[C] lambdaPrime = exp(log_lambdaPrime);
    vector<lower=0>[C] NSERate     = exp(log_NSERate);
    // Noncentered log_phi with a sum-to-zero z_phi, GEOMETRIC-mean-one variant.
    // phi_use_data=0: mean(log_phi) = 0 EXACTLY (geometric mean phi = 1), so the
    //   phi LOCATION is independent of sphi (sphi is pure spread).  This breaks
    //   the sphi<->global-lambda ridge that the -0.5*sphi^2 offset induces in
    //   panse_sphi_est_sumzero.stan (raising sphi there lowers mean(log_phi),
    //   which the global lambda scale must track -> cor(sphi,lambda) ~ 0.80).
    // phi_use_data=1: per-gene prior mean/SD; deviations sum to zero.
    vector[G] log_phi;
    if (phi_use_data) {
        log_phi = phi_prior_mu + phi_prior_sigma .* z_phi;
    } else {
        log_phi = sphi * z_phi;   // geometric-mean-one (median phi=1); no offset
    }
    vector<lower=0>[G] phi = exp(log_phi);

    vector[C] log_alpha_term;
    vector[C] log_psuccess;
    for (c in 1:C) {
        log_alpha_term[c] = log_alpha[c] - log_U - log_lambdaPrime[c];
        log_psuccess[c] = log_psuccess_hybrid(alpha[c], lambdaPrime[c], NSERate[c]);
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

    /* NO soft mean(phi)=1 anchor: the sum_to_zero_vector pins mean(z_phi)=0 by
     * construction, so the phi scale gauge is removed exactly (the hard
     * reparameterization analog of native fix.z). */

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
