/* ============================================================================
 * roc_mixture_sphi.stan -- ROC HMC with 2-component lognormal phi prior.
 *
 * Mirrors the phi-mixture-LN prior implemented in C++ on
 * feat/rst-mixture-state (task #12 family).  Per gene:
 *
 *   log_phi_g | mixture-LN  ~  p * Normal(mu1, sigma1) + (1 - p) * Normal(mu2, sigma2)
 *
 * with mu2 DERIVED from the constraint mean(phi) = 1 (PHI_CONSTRAINT_MEAN):
 *
 *   p * exp(mu1 + sigma1^2/2) + (1 - p) * exp(mu2 + sigma2^2/2) = 1
 *   => mu2 = log( (1 - p * exp(mu1 + sigma1^2/2)) / (1 - p) ) - sigma2^2/2
 *
 * Label-switching guard: require mu2 >= mu1 (component 2 has the higher
 * mean in log-space).  Together with mean=1 constraint this also requires
 * the closed-form numerator (1 - p * exp(mu1 + sigma1^2/2)) > 0; otherwise
 * the constraint is infeasible and the proposal is rejected.
 *
 * Free hyperparameters: p, mu1, sigma1, sigma2.  mu2 is a transformed
 * parameter.  Hyperpriors match the C++ defaults in src/Parameter.cpp:
 *
 *   p      ~ Beta(8, 2)            (mean 0.8; concentrated toward larger
 *                                    weight on component 1)
 *   mu1    ~ Normal(0, 10)         (very weak; matches phi mixture hyper
 *                                    Mu1_mean=0, Mu1_sd=10)
 *   sigma1 ~ half_normal(0, 1)     (matches sigma1_scale=1 default)
 *   sigma2 ~ half_normal(0, 1)     (matches sigma2_scale=1 default)
 *
 * The half_normal vs half_cauchy choice for the sigma_k priors is a
 * judgment call; both have scale=1 as the C++ default.  half_normal is
 * lighter-tailed and tends to mix better in HMC; switch to cauchy via
 * sigma ~ cauchy(0, 1) if posterior tails matter.
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
    vector[G] log_phi;
    real<lower=0, upper=1> p;
    real mu1;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters {
    vector<lower=0>[G] phi = exp(log_phi);

    // Derived mu2 from the mean=1 constraint.
    // Infeasibility (numerator <= 0) is handled in the model block by
    // adding -infinity to the log-posterior; here we just compute it.
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
    // If numer <= 0, mu2 is not well-defined; reject.  Likewise for mu2 < mu1
    // (component 2 must dominate in log-space).
    if (numer <= 0)  reject("mixture-LN constraint infeasible: numer <= 0");
    if (mu2 < mu1)   reject("label-switching guard: mu2 < mu1");

    // Per-gene mixture-LN log prior on log_phi
    for (g in 1:G) {
        target += log_mix(
            p,
            normal_lpdf(log_phi[g] | mu1, sigma1),
            normal_lpdf(log_phi[g] | mu2, sigma2));
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
