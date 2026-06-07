/* ============================================================================
 * panse_sphi_est_noncentered_invgamma.stan -- INVERSE-GAMMA dwell variant of
 * panse_sphi_est_noncentered.stan.
 *
 * DWELL DISTRIBUTION: the per-codon dwell time W ~ InverseGamma(shape=alpha,
 * scale=lambda).  Equivalently the per-codon elongation RATE 1/W ~
 * Gamma(shape=alpha, rate=lambda) -- a Gamma on the rate rather than the
 * waiting time.  Two consequences vs the Gamma-dwell parent:
 *
 *   (1) SURVIVAL  -- EXACT, reuses the parent's incomplete-gamma machinery.
 *       Substituting u = 1/w turns W ~ InvGamma(a,l) into U=1/W ~ Gamma(a,
 *       rate=l), and
 *         psuccess = E_W[v/(W+v)]
 *                  = 1 - E_{U~Gamma(a,rate=l)}[ NSERate/(U+NSERate) ].
 *       The inner expectation is the SAME functional the Gamma model uses for
 *       its survival, only with the role of "v" played by vv=NSERate (instead
 *       of v=1/NSERate).  So we reuse log_upper_incomplete_gamma verbatim:
 *         inner = exp( a*log(l*vv) + l*vv + logGammaUpper(1-a, l*vv) ), vv=NSERate
 *         log_psuccess = log1m(inner).
 *       Bounded in (0,1), monotone decreasing in NSERate; ->1 as NSERate->0.
 *       Verified to ~1e-16 vs direct integration.
 *
 *   (2) COUNTS  -- the marginal is the SICHEL distribution (Poisson-InvGamma),
 *       whose PMF involves K_nu(.) (modified Bessel of the 2nd kind) of order
 *       nu=y-alpha.  Stan 2.x has no Bessel-K, and K_nu OVERFLOWS for the large
 *       orders that real footprint counts (hundreds-thousands) imply, so a
 *       faithful Sichel is numerically fragile and silently wrong at data
 *       scale.  We therefore ship a DOCUMENTED MOMENT-MATCHED SURROGATE.
 *
 *       The Sichel mean and variance are EXACTLY
 *         E[Y]   = mu
 *         Var[Y] = mu + mu^2/(alpha-2)        (requires alpha>2)
 *       i.e. the NB2 mean-variance form with size = alpha-2 (vs the Gamma
 *       model's size = alpha).  We match those two moments:
 *         Y ~ NB2(mean=mu, size = alpha-2).
 *       This is finite everywhere, free, identical in form to the Gamma count
 *       model, and exact in mean+variance.  It is an APPROXIMATION: the true
 *       Sichel has heavier (Bessel-K) tails, so the surrogate slightly
 *       under-weights extreme counts.  For LOO comparison of dwell-distribution
 *       ASSUMPTIONS this is the right call (the distinctive survival term is
 *       exact; the count term differs from Gamma only in dispersion, which is
 *       exactly the moment-level signal LOO resolves).  See
 *       notes/dwell-time-distributions.md section 3 for the derivation,
 *       the overflow evidence, and the faithful-vs-surrogate recommendation.
 *
 * GUARD: size = alpha-2 requires alpha>2 for a finite Sichel variance.  When
 * alpha<=2 the InvGamma dwell has no finite variance and the surrogate is
 * undefined; we clamp the NB size to a small positive floor (SIZE_FLOOR) so
 * the sampler does not hit a domain error while exploring, and rely on the
 * data to keep alpha>2 in the bulk.  (The log_alpha lower bound can also be
 * raised to log(2) via config if a hard floor is preferred.)
 *
 * The parameter vector is byte-identical to panse_sphi_est_noncentered.stan
 * (log_alpha[C], log_lambdaPrime[C], log_NSERate[C], sphi, z_phi[G]).
 * Everything else (noncentered log_phi, generalized phi prior, soft mean(phi)=1
 * anchor, partial_sum, generated quantities) is identical to the parent.
 * ============================================================================ */

functions {
    /* log of the UNregularized upper incomplete gamma Gamma(s, x), x > 0, via
     * the forward modified-Lentz continued fraction (identical to the parent
     * Gamma model; see panse_sphi_est_noncentered.stan header for the math and
     * the q* breach history).  Stable for all real s including s<=0. */
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

    /* InverseGamma-dwell log survival:
     *   psuccess = 1 - E_{U~Gamma(alpha,rate=lambda)}[ NSERate/(U+NSERate) ].
     * The inner expectation reuses the Gamma-survival closed form with the
     * "v" argument set to vv=NSERate (NOT 1/NSERate).  log_psuccess = log1m(inner).
     * For very small NSERate, inner is tiny and log1m(inner) is well
     * conditioned; we cap inner just below 1 to keep log1m finite.
     *
     * NUMERICS: the inner expectation's incomplete-gamma argument is
     * lvv = lambda*NSERate, which is SMALL in the operating regime (lambda in
     * [1e-3,100], NSERate in [1e-7,1e-3] default -> lvv <= 0.1; wide-NSE
     * NSERate up to 0.1 -> lvv up to 10).  The forward modified-Lentz CF (used
     * by the parent Gamma model) converges fast only for LARGE argument; here
     * the argument is small, the OPPOSITE corner, where the CF is slow / fails
     * to converge in 2000 iters.  So we split by lvv:
     *   lvv < CF_SWITCH (0.5): 2-term asymptotic of inner = E_U[vv/(U+vv)]
     *       inner ~ vv*E[1/U] - vv^2*E[1/U^2], U~Gamma(a,rate=l), vv=NSERate
     *             = lvv/(a-1) - lvv^2/((a-1)(a-2)).
     *       Verified: max |log-survival error| ~ 5e-6 over the full default box
     *       (a in [2,100], lvv in [1e-10,0.1]) vs direct integration.  Requires
     *       a>1 (E[1/U] finite); for a<=1 the asymptotic breaks but survival
     *       there is ~1 anyway (inner tiny), so we floor a-1 below.
     *   lvv >= CF_SWITCH: the CF converges fast (worst ~70 iters), use it.
     */
    real log_psuccess_invgamma(real alpha, real lambda, real NSERate) {
        real CF_SWITCH = 0.5;
        real vv  = NSERate;
        real lvv = lambda * vv;
        real log_inner;
        if (lvv < CF_SWITCH) {
            // 2-term small-argument asymptotic of the inner expectation.
            real am1 = fmax(alpha - 1, 1e-6);       // guard E[1/U] for alpha<=1
            real am2 = fmax(alpha - 2, 1e-6);       // guard E[1/U^2] for alpha<=2
            real inner_lin = lvv / am1 - lvv * lvv / (am1 * am2);
            log_inner = log(fmax(inner_lin, 1e-300));
        } else {
            log_inner = alpha * log(lvv) + lvv + log_upper_incomplete_gamma(1 - alpha, lvv);
        }
        real inner = exp(log_inner);
        if (inner >= 1) inner = 1 - 1e-12;          // guard log1m at the boundary
        return log1m(inner);
    }

    real partial_sum(array[] int slice_g, int start, int end,
                     array[] int gene_offset,
                     array[] int codon_at_pos,
                     array[] int y,
                     array[] int like_mask,
                     int all_unmasked,
                     vector nb_size,                 // moment-matched NB2 size = alpha-2 (floored)
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
                vector[n] size_g;
                for (j in 1:n) {
                    int c = codon_at_pos[p0 + j - 1];
                    log_psuccess_g[j] = log_psuccess[c];
                    size_g[j]         = nb_size[c];
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
                lp += neg_binomial_2_log_lpmf(y[p0:p1] | log_mu_g, size_g);
            } else {
                real log_survive = 0;
                for (p in p0:p1) {
                    int c = codon_at_pos[p];
                    if (like_mask[p] == 1) {
                        lp += neg_binomial_2_log_lpmf(
                            y[p] | log_alpha_term[c] + lpg + log_survive,
                            nb_size[c]);
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
    real SIZE_FLOOR = 1e-2;   // floor for NB size = alpha-2 when alpha<=2
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
    vector<lower=0>[C] nb_size;        // moment-matched Sichel surrogate size
    for (c in 1:C) {
        log_alpha_term[c] = log_alpha[c] - log_U - log_lambdaPrime[c];
        // Exact InvGamma survival (reuses the Gamma incomplete-gamma machinery).
        log_psuccess[c] = log_psuccess_invgamma(alpha[c], lambdaPrime[c], NSERate[c]);
        // Sichel-count surrogate: NB2 size = alpha-2, floored for alpha<=2.
        nb_size[c] = fmax(alpha[c] - 2, SIZE_FLOOR);
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
                         nb_size, log_alpha_term, log_psuccess, log_phi);
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
                        nb_size[c]);
                } else {
                    log_lik[p] = 0;
                }
                log_survive += log_psuccess[c];
            }
        }
    }
}
