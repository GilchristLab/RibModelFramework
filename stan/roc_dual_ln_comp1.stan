/* ============================================================================
 * roc_dual_ln_comp1.stan
 *
 * Component-1-anchored mixture-LN ROC model (task #14).
 *
 * The unit-median(phi) constraint is applied to the BULK component ONLY,
 * not to the full mixture:
 *
 *     mphi1 = -sphi1^2 / 2      (=> median(phi | comp 1) = 1)
 *
 * This makes the bulk's (mphi1, sphi1) directly comparable to a single-LN
 * fit on the same data, and lets a single-LN posterior on sphi serve as
 * an informative prior on sphi1 (the bulk SD) under matched-data conditions.
 *
 * The tail component (component 2) is anchored relative to a single-LN
 * posterior summary, passed in via data:
 *
 *     q_log = mphi_s + Phi^{-1}(tail_quantile) * sphi_s
 *
 * where (mphi_s, sphi_s) = posterior mean of (mphi, sphi) from a single-LN
 * fit on the same data.  Hard bounds:
 *
 *     mphi2  >= q_log
 *     sphi2  <= sphi1 / ratio_min   (sampled via sphi_ratio in (0, 1/ratio_min])
 *     p      >= p_min               (>= 1 - p is the tail weight)
 *
 * See analysis/02x_Dual.LN.ROC.MCMC.Fits/notes/component1-anchored-dual-LN.md
 * (Lokiarchaeota repo) for design rationale and YAML schema.
 *
 * Naming convention: mphi / sphi throughout (task #15 phase 1, 2026-05-28).
 * Backends: Stan (this file) and (native) MCMC.  AnaCoDa is the C++ library
 * underneath native MCMC, not a separate backend.
 *
 * Centered / non-centered toggle via `noncentered` data flag, both anchored
 * on component 1:
 *
 *   noncentered = 0 (default):
 *     latent_phi[g] IS log_phi[g].  Mixture prior in log_phi-space.
 *
 *   noncentered = 1:
 *     latent_phi[g] is z_phi anchored on component 1.
 *     log_phi[g] = mphi1 + sphi1 * latent_phi[g]    (transformed parameter)
 *     Mixture prior in z-space:
 *       z_phi ~ p * N(0, 1) + (1 - p) * N(delta, ratio)
 *     where delta = (mphi2 - mphi1) / sphi1,  ratio = sphi2 / sphi1.
 *     Removes the (log_phi, sphi1) funnel for the bulk; tail funnel
 *     partially mitigated.
 *
 * THREADING: compile with cpp_options = list(stan_threads = TRUE) and pass
 * threads_per_chain > 1 to enable reduce_sum parallelism.
 * ============================================================================ */

functions {
    real partial_sum(array[] int slice_g, int start, int end,
                          int A,
                          array[] int aa_start, array[] int aa_end,
                          array[,] int y_k,
                          array[,] int N_ga,
                          vector dM, vector dEta,
                          vector phi) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);
            for (a in 1:A) {
                if (N_ga[g, a] == 0) continue;
                int s = aa_start[a];
                int e = aa_end[a];
                int n = e - s + 1;
                vector[n + 1] eta_full;
                eta_full[1] = 0;
                for (k in 1:n) {
                    eta_full[k + 1] = -dM[s - 1 + k] - dEta[s - 1 + k] * phi[g];
                }
                lp += -N_ga[g, a] * log_sum_exp(eta_full);
            }
        }
        return lp;
    }

    // Centered mixture prior on log_phi.
    real partial_mix_prior(array[] int slice_g, int start, int end,
                                vector log_phi,
                                real p, real mphi1, real mphi2,
                                real sphi1, real sphi2) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += log_mix(
                p,
                normal_lpdf(log_phi[g] | mphi1, sphi1),
                normal_lpdf(log_phi[g] | mphi2, sphi2));
        }
        return lp;
    }

    // Non-centered mixture prior on z_phi (component 1 anchored at N(0, 1)).
    real partial_mix_prior_z(array[] int slice_g, int start, int end,
                                  vector z_phi,
                                  real p, real delta, real ratio) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            lp += log_mix(
                p,
                normal_lpdf(z_phi[g] | 0, 1),
                normal_lpdf(z_phi[g] | delta, ratio));
        }
        return lp;
    }
}

data {
    int<lower=1> G;
    int<lower=1> A;
    int<lower=1> K;

    array[A] int<lower=1> aa_start;
    array[A] int<lower=1> aa_end;

    array[G, K] int<lower=0> y_k;
    array[G, A] int<lower=0> N_ga;

    vector[K] dM_prior_mean;
    vector<lower=0>[K] dM_prior_sd;
    vector[K] dEta_prior_mean;
    vector<lower=0>[K] dEta_prior_sd;

    // ---- Single-LN posterior reference (wrapper precomputes from manifest) ----
    real q_log;                            // mphi_s + qnorm(tail_quantile) * sphi_s
    real<lower=0> mphi2_prior_scale;       // half-normal scale on (mphi2 - q_log)

    // ---- Bulk-component (component 1) hyperprior ----
    real<lower=0> sphi1_prior_scale;       // half-normal scale on sphi1

    // ---- Scale-ratio constraint: sphi_ratio = sphi2 / sphi1 ----
    real<lower=0, upper=1> inv_ratio_min;  // = 1 / ratio_min (e.g. 0.2 for ratio_min=5)
    real<lower=0> sphi_ratio_prior_scale;  // half-normal scale on sphi_ratio

    // ---- Bulk weight constraint: p in [p_min, 1] ----
    real<lower=0, upper=1> p_min;          // e.g. 0.95
    real<lower=0> p_alpha;                 // Beta shape on p (truncated to [p_min, 1])
    real<lower=0> p_beta;

    int<lower=0, upper=1> noncentered;     // 0 = centered (default), 1 = non-centered
    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    // latent_phi[g] = log_phi[g] when centered, z_phi[g] when non-centered.
    vector[G] latent_phi;
    real<lower=p_min, upper=1.0> p;
    real<lower=0> sphi1;
    real<lower=0, upper=inv_ratio_min> sphi_ratio;
    real<lower=q_log> mphi2;
}

transformed parameters {
    // Component-1 anchor: median(phi | comp 1) = 1.
    real mphi1 = -0.5 * sphi1 * sphi1;
    real sphi2 = sphi1 * sphi_ratio;

    // log_phi: identity when centered; derived from z_phi when non-centered.
    vector[G] log_phi = noncentered ? (mphi1 + sphi1 * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    // z-space mixture parameters (used only when noncentered=1; computed
    // unconditionally to avoid branching in transformed parameters).
    real delta = (mphi2 - mphi1) / sphi1;
    real ratio = sphi_ratio;
}

model {
    // Hyperpriors.  Bounded declarations produce truncated densities; the
    // truncation factors are constants and do not affect MCMC posterior shape.
    p          ~ beta(p_alpha, p_beta);
    sphi1      ~ normal(0, sphi1_prior_scale);          // half-normal via <lower=0>
    sphi_ratio ~ normal(0, sphi_ratio_prior_scale);     // half-normal on (0, inv_ratio_min]
    mphi2      ~ normal(q_log, mphi2_prior_scale);      // truncated to mphi2 >= q_log

    // CSP priors (vectorized; outside reduce_sum since vector ops already
    // benefit from Stan's vectorization).
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Mixture prior on per-gene log_phi (threaded).
    if (noncentered) {
        target += reduce_sum(partial_mix_prior_z, gene_indices, grainsize,
                             latent_phi, p, delta, ratio);
    } else {
        target += reduce_sum(partial_mix_prior, gene_indices, grainsize,
                             latent_phi, p, mphi1, mphi2, sphi1, sphi2);
    }

    // Per-gene multinomial-logit likelihood (threaded).
    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta, phi);
}
