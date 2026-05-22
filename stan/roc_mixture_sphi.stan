/* ============================================================================
 * roc_mixture_sphi.stan -- ROC HMC with 2-component lognormal phi prior.
 *
 * Non-centered on component 1 (the dominant component under the default
 * p ~ Beta(8, 2) prior, which puts E[p] = 0.8 on component 1).
 *
 *   z_phi[g] ~ std_normal()                (parameter)
 *   log_phi[g] = mu1 + sigma1 * z_phi[g]   (transformed parameter)
 *
 * Under this change of variables, the original mixture-LN prior on log_phi
 *
 *   log_mix(p, N(log_phi | mu1, sigma1), N(log_phi | mu2, sigma2))
 *
 * becomes, in z_phi-space (with Jacobian |d log_phi / d z_phi| = sigma1):
 *
 *   log_mix(p,
 *           std_normal_lpdf(z_phi),
 *           N(log_phi | mu2, sigma2) + log(sigma1))
 *
 * Component 1's contribution is exactly std_normal (no funnel).  Component 2's
 * contribution is still centered but its weight is small (1 - E[p] ~ 0.2),
 * so the residual funnel is rare in posterior draws.  For workloads where
 * component 2 dominates (atypical), the analogous "non-centered on
 * component 2" variant would be the right tool; not implemented here.
 *
 * Mirrors the phi-mixture-LN prior implemented in C++ on feat/rst-mixture-state
 * (task #12 family).  mu2 is DERIVED from mean(phi) = 1 (PHI_CONSTRAINT_MEAN):
 *
 *   p * exp(mu1 + sigma1^2/2) + (1 - p) * exp(mu2 + sigma2^2/2) = 1
 *   => mu2 = log( (1 - p * exp(mu1 + sigma1^2/2)) / (1 - p) ) - sigma2^2/2
 *
 * Label-switching guard: mu2 >= mu1 (component 2 is the higher-mode
 * lognormal in log-space).  Constraint infeasibility (numer <= 0) reject.
 *
 * Hyperpriors match C++ defaults in src/Parameter.cpp:
 *   p      ~ Beta(p_alpha, p_beta)         default Beta(8, 2)
 *   mu1    ~ Normal(mu1_prior_mean, mu1_prior_sd)   default N(0, 10)
 *   sigma1 ~ half_normal(0, sigma1_prior_scale)     default scale 1
 *   sigma2 ~ half_normal(0, sigma2_prior_scale)     default scale 1
 *
 * See roc_basic.stan for the data layout (identical here; only the prior
 * on log_phi changes).
 * ============================================================================ */

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

    // Hyperprior knobs (override the defaults set in this file's docstring).
    real<lower=0> p_alpha;
    real<lower=0> p_beta;
    real         mu1_prior_mean;
    real<lower=0> mu1_prior_sd;
    real<lower=0> sigma1_prior_scale;
    real<lower=0> sigma2_prior_scale;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[G] z_phi;                        // non-centered latent for component 1
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    // log_phi is deterministic in (mu1, sigma1, z_phi); funnel-free in
    // sampling space for component-1 (the dominant component).
    vector[G] log_phi = mu1 + sigma1 * z_phi;
    vector<lower=0>[G] phi = exp(log_phi);

    // Derived mu2 from the mean=1 constraint.
    // Infeasibility (numerator <= 0) is handled in the model block by
    // rejecting; here we just compute it (Stan needs no try/catch).
    real numer = 1.0 - p * exp(mu1 + 0.5 * sigma1 * sigma1);
    real mu2 = log(numer / (1.0 - p)) - 0.5 * sigma2 * sigma2;
}

model {
    // Hyperpriors
    p      ~ beta(p_alpha, p_beta);
    mu1    ~ normal(mu1_prior_mean, mu1_prior_sd);
    sigma1 ~ normal(0, sigma1_prior_scale);             // half-normal via lower=0
    sigma2 ~ normal(0, sigma2_prior_scale);

    // Constraint feasibility + label-switching guard
    if (numer <= 0) reject("mixture-LN constraint infeasible: numer <= 0");
    if (mu2 < mu1)  reject("label-switching guard: mu2 < mu1");

    // Per-gene mixture-LN log prior on log_phi, expressed in z_phi-space.
    // Change of variables (log_phi = mu1 + sigma1 * z_phi):
    //   component 1: N(log_phi | mu1, sigma1) + log|sigma1| = std_normal_lpdf(z_phi)
    //   component 2: N(log_phi | mu2, sigma2) + log|sigma1|
    // The log(sigma1) Jacobian appears only on the component-2 leg because
    // it cancels exactly on the component-1 leg.
    for (g in 1:G) {
        target += log_mix(
            p,
            std_normal_lpdf(z_phi[g]),
            normal_lpdf(log_phi[g] | mu2, sigma2) + log(sigma1));
    }

    // dM, dEta priors (Gaussian, per-codon means/SDs)
    dM   ~ normal(dM_prior_mean, dM_prior_sd);
    dEta ~ normal(dEta_prior_mean, dEta_prior_sd);

    // Likelihood (identical to roc_basic.stan)
    for (g in 1:G) {
        target += dot_product(to_vector(y_k[g, :]), -dM - dEta * phi[g]);
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
            target += -N_ga[g, a] * log_sum_exp(eta_full);
        }
    }
}
