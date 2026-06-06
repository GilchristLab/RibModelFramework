/* ============================================================================
 * roc_sphi_est.stan -- ROC HMC, sphi estimated, centered + non-centered,
 *                      reduce_sum threaded.  phi.mphi / phi.sphi spec API.
 *
 * Parameterisation switch (data field `noncentered`):
 *   noncentered = 0  (centered, default):
 *     latent_phi[g] IS log_phi[g].
 *     Prior: log_phi ~ Normal(mphi, sphi).
 *     Best when data strongly anchors each gene (G >= 1000, full-genome).
 *
 *   noncentered = 1:
 *     latent_phi[g] is z_phi[g] ~ N(0,1); log_phi = mphi + sphi * z_phi.
 *     Removes Neal's funnel for data-sparse fits with poor sphi ESS.
 *
 * phi.mphi spec (data fields phi_mphi_mode / phi_mphi_statistic / phi_mphi_value
 *                             / phi_mphi_fixed / phi_mphi_prior_type /
 *                             phi_mphi_prior_mean / phi_mphi_prior_sd):
 *   phi_mphi_mode = 0  (constrained):
 *     mphi is derived from sphi each iteration so that the chosen statistic
 *     of phi equals phi_mphi_value.  Statistics (phi_mphi_statistic):
 *       0 = mean:     mphi = log(v) - sphi^2/2
 *       1 = median:   mphi = log(v)
 *       2 = mode:     mphi = log(v) + sphi^2
 *       3 = variance: mphi = (log(v / (exp(sphi^2) - 1)) - sphi^2) / 2
 *       4 = sd:       mphi = (log(v^2 / (exp(sphi^2) - 1)) - sphi^2) / 2
 *
 *   phi_mphi_mode = 1  (fixed):
 *     mphi = phi_mphi_fixed  (a data constant; sphi is still sampled).
 *
 *   phi_mphi_mode = 2  (estimated):
 *     mphi_param is a free sampled parameter.  Prior controlled by
 *     phi_mphi_prior_type (0=flat, 1=Normal(mean,sd)).
 *     Enables deta_scale_anchor to actually collapse the dEta-phi ridge.
 *
 * phi.sphi spec (data fields sphi_low / sphi_high / sphi_prior_type /
 *                             sphi_prior_mean / sphi_prior_sd):
 *   sphi is declared real<lower=sphi_low, upper=sphi_high>.
 *   sphi_prior_type = 0: bounds-only (Uniform implied; no prior statement).
 *   sphi_prior_type = 1: sphi ~ normal(sphi_prior_mean, sphi_prior_sd)
 *                        truncated to positive by the lower=0 bound.
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
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

    // phi.sphi spec
    real<lower=0> sphi_low;               // lower bound for sphi
    real          sphi_high;              // upper bound (large for half-line)
    int<lower=0, upper=1> sphi_prior_type; // 0=bounds-only (uniform), 1=normal
    real          sphi_prior_mean;        // used when sphi_prior_type=1
    real<lower=0> sphi_prior_sd;          // used when sphi_prior_type=1

    // phi.mphi spec
    int<lower=0, upper=2> phi_mphi_mode;      // 0=constrained, 1=fixed, 2=estimated
    int<lower=0, upper=4> phi_mphi_statistic; // 0=mean,1=median,2=mode,3=variance,4=sd
    real<lower=0>         phi_mphi_value;     // target statistic value (constrained mode)
    real                  phi_mphi_fixed;     // fixed mphi value (fixed mode)
    // phi.mphi prior (used when phi_mphi_mode==2; placeholder values for modes 0/1)
    int<lower=0, upper=1> phi_mphi_prior_type; // 0=flat, 1=Normal
    real                  phi_mphi_prior_mean;
    real<lower=0>         phi_mphi_prior_sd;

    int<lower=0, upper=1> noncentered;        // 0=centered (default), 1=non-centered

    // dEta / phi decorrelation flags (orthogonal to phi spec)
    int<lower=0, upper=1> deta_scale_anchor;  // 0=off (default), 1=dS reparameterisation
    real                  deta_anchor_ref;    // reference mphi for scale anchor
    real                  deta_phi_center;    // phi centering constant (0=off)

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
    real<lower=sphi_low, upper=sphi_high> sphi;
    // mphi_param: free parameter when phi_mphi_mode==2; pinned via std_normal() otherwise.
    real mphi_param;
}

transformed parameters {
    // Compute mphi from phi spec
    real s2   = square(sphi);
    real mphi;
    if (phi_mphi_mode == 2) {
        mphi = mphi_param;                                 // estimated (free parameter)
    } else if (phi_mphi_mode == 1) {
        mphi = phi_mphi_fixed;                             // fixed constant
    } else if (phi_mphi_statistic == 0) {
        mphi = log(phi_mphi_value) - 0.5 * s2;            // constrained: mean
    } else if (phi_mphi_statistic == 1) {
        mphi = log(phi_mphi_value);                        // constrained: median
    } else if (phi_mphi_statistic == 2) {
        mphi = log(phi_mphi_value) + s2;                   // constrained: mode
    } else if (phi_mphi_statistic == 3) {
        mphi = 0.5 * (log(phi_mphi_value / (exp(s2) - 1.0)) - s2);  // constrained: variance
    } else {
        mphi = 0.5 * (log(square(phi_mphi_value) / (exp(s2) - 1.0)) - s2);  // constrained: sd
    }

    // log_phi: identity when centered; derived from z_phi when non-centered.
    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    // Effective selection (scale-anchor reparameterisation).
    vector[K] dEta_eff = (deta_scale_anchor == 1)
        ? dEta * exp(-(mphi - deta_anchor_ref))
        : dEta;

    // Phi centering: decorrelates dM and dEta by shifting the predictor.
    vector[G] phi_eff = phi - deta_phi_center;
    vector[K] dM_at0  = dM - dEta_eff * deta_phi_center;
}

model {
    // Prior on the intercept-at-phi=0 (dM_at0 == dM when deta_phi_center=0).
    dM_at0 ~ normal(dM_prior_mean, dM_prior_sd);
    dEta   ~ normal(dEta_prior_mean, dEta_prior_sd);

    // sphi prior: bounds always enforced by the parameter constraint.
    // A normal prior statement is added only when sphi_prior_type=1.
    if (sphi_prior_type == 1) {
        sphi ~ normal(sphi_prior_mean, sphi_prior_sd);
    }

    // mphi_param: prior when estimated (mode==2); pinned via std_normal() otherwise.
    // When mode != 2, mphi_param is independent of all likelihood terms, so
    // std_normal() is a harmless regulariser with negligible HMC overhead.
    if (phi_mphi_mode == 2) {
        if (phi_mphi_prior_type == 1) {
            mphi_param ~ normal(phi_mphi_prior_mean, phi_mphi_prior_sd);
        }
    } else {
        mphi_param ~ std_normal();
    }

    // Phi prior
    if (noncentered) {
        latent_phi ~ std_normal();
    } else {
        latent_phi ~ normal(mphi, sphi);
    }

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, y_k, N_ga,
                         dM, dEta_eff, phi_eff);
}
