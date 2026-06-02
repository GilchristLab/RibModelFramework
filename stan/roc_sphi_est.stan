/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC, sphi estimated, centered + non-centered,
 *                      reduce_sum threaded.
 *
 * Unified file replacing roc_sphi_est_centered.stan and
 * roc_sphi_est_noncentered.stan.  The data field `noncentered` selects
 * the parameterization at run time without recompilation:
 *
 *   noncentered = 0  (centered, default):
 *     latent_phi[g] IS log_phi[g].
 *     Prior: log_phi ~ Normal(mphi, sphi)  where mphi = -sphi^2/2
 *     Best when data strongly anchors each gene (G >= 1000, full-genome fits).
 *
 *   noncentered = 1:
 *     latent_phi[g] is z_phi[g] ~ N(0,1).
 *     log_phi[g] = mphi + sphi * z_phi[g]   (transformed parameter)
 *     Removes Neal's funnel in (log_phi, sphi) for data-sparse fits.
 *     Use when centered has poor sphi ESS or chains stuck at small sphi.
 *
 * Phi-scale anchor: data flag `anchor_phi` selects the convention.
 *   anchor_phi = 0  (default, backward compatible):
 *     mphi = -sphi^2/2  (mean(phi) = 1 in expectation, v.3 convention)
 *   anchor_phi = 1:
 *     mphi is a sampled parameter (mphi_param) with soft prior
 *       mphi ~ Normal(0, mphi_prior_sd)
 *     so median(phi) = exp(mphi) is softly anchored near 1.  Recommended
 *     for cross-genome phi comparability; also collapses the (dEta, phi)
 *     scale ridge that mean-anchored sphi alone leaves under-determined.
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
 * threads_per_chain > 1 to $sample() to enable reduce_sum parallelism.
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
    real          sphi_prior_mean;      // prior mean for sphi (default 0 = half-normal at 0)
    real<lower=0> sphi_prior_sd;        // SD for sphi prior; lower=0 truncates to positive

    int<lower=0, upper=1> noncentered;  // 0=centered (default), 1=non-centered
    int<lower=0, upper=1> anchor_phi;   // 0=mean(phi)=1 via mphi=-sphi^2/2 (default), 1=median(phi)=1 via soft prior on mphi
    real<lower=0> mphi_prior_sd;        // SD of soft median anchor prior; ignored when anchor_phi=0

    // Scale-anchored selection (phi-dEta ridge reparameterization).
    //   deta_scale_anchor = 0 (default): dEta enters the likelihood directly.
    //   deta_scale_anchor = 1: the sampled `dEta` vector is interpreted as dS
    //     (selection at the anchor phi level) and the likelihood uses
    //     dEta_eff = dS * exp(-(mphi - deta_anchor_ref)).  The likelihood only
    //     ever uses the product dEta*phi and phi = exp(mphi + dev), so the
    //     global level mphi cancels out of that product and dS is decorrelated
    //     from mphi -- collapsing the dominant dEta-phi scale ridge.  REQUIRES
    //     anchor_phi=1 so mphi is an explicit (sampled) parameter tracking the
    //     global phi level (with anchor_phi=0, mphi=-sphi^2/2 does not track
    //     the free mean(log_phi) drift, so the ridge is only partly removed).
    int<lower=0, upper=1> deta_scale_anchor;
    real deta_anchor_ref;               // reference mphi level dS anchors at; ignored when deta_scale_anchor=0

    // Phi-centering (dM-dEta intercept-slope decorrelation).  The likelihood
    //   eta_{g,c} = -dM_c - dEta_c * phi_g
    // is a per-codon regression of eta on phi with intercept -dM_c and slope
    // -dEta_c; an uncentered predictor (phi_g ~ 1-4, not 0) makes intercept
    // and slope strongly collinear (measured corr(dM_c,dEta_c) ~ -0.4 on every
    // codon).  Subtracting a center c decorrelates them:
    //   eta_{g,c} = -dM_c - dEta_eff_c * (phi_g - c)
    // Here the sampled `dM` is the intercept AT phi=c; the dM prior is applied
    // to the reconstructed intercept-at-0, dM_at0 = dM - dEta_eff*c, so an
    // informative dM prior (scuo_low/encp_low) keeps its meaning and is not
    // biased by c.  deta_phi_center = 0 disables centering (backward compatible).
    // Composes independently with deta_scale_anchor.
    real deta_phi_center;
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
    real<lower=0> sphi;
    real mphi_param;                    // used only when anchor_phi=1; phantom otherwise
}

transformed parameters {
    // mphi: derived from sphi (anchor_phi=0, mean=1 convention) or sampled
    // via mphi_param (anchor_phi=1, soft median=1 anchor).
    real mphi = (anchor_phi == 1) ? mphi_param : -0.5 * sphi * sphi;
    // log_phi: identity when centered; derived from z_phi when non-centered.
    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    // Effective selection used in the likelihood.  When deta_scale_anchor=1
    // the sampled `dEta` is dS and dEta_eff rescales it by exp(-(mphi-ref))
    // so the mphi dependence cancels in the dEta_eff*phi product (see data
    // block).  When 0, dEta_eff == dEta (backward compatible).
    vector[K] dEta_eff = (deta_scale_anchor == 1)
        ? dEta * exp(-(mphi - deta_anchor_ref))
        : dEta;

    // Phi-centering: predictor used in the likelihood, and the intercept-at-0
    // (dM_at0) that the dM prior is placed on.  deta_phi_center=0 => phi_eff
    // == phi and dM_at0 == dM (identical to the un-centered model).
    vector[G] phi_eff = phi - deta_phi_center;
    vector[K] dM_at0  = dM - dEta_eff * deta_phi_center;
}

model {
    // Prior on the intercept-at-phi=0 so an informative dM prior keeps its
    // meaning under centering (dM_at0 == dM when deta_phi_center=0).  The map
    // (dM,dEta) -> (dM_at0,dEta) is a unit-Jacobian shear, so no correction.
    dM_at0 ~ normal(dM_prior_mean, dM_prior_sd);
    dEta   ~ normal(dEta_prior_mean, dEta_prior_sd);
    sphi ~ normal(sphi_prior_mean, sphi_prior_sd);  // truncated to positive via lower=0 constraint

    // Phi-scale anchor prior on mphi_param.  When anchor_phi=1, this softly
    // pins median(phi) = exp(mphi) near 1.  When anchor_phi=0, mphi_param
    // is a phantom (unused in transformed parameters) -- we still give it a
    // wide proper prior so the posterior stays proper.
    if (anchor_phi == 1) {
        mphi_param ~ normal(0, mphi_prior_sd);
    } else {
        mphi_param ~ normal(0, 1);      // phantom; never enters mphi or likelihood
    }

    // Phi prior: centered form uses the lognormal prior directly on log_phi;
    // non-centered form places std_normal on the latent z_phi.
    if (noncentered) {
        latent_phi ~ std_normal();
    } else {
        latent_phi ~ normal(mphi, sphi);
    }

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta_eff, phi_eff);
}
