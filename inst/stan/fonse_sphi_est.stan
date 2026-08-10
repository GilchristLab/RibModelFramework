/* ============================================================================
 * fonse_sphi_est.stan -- FONSE HMC, sphi estimated, centered + non-centered,
 *                        reduce_sum threaded.  phi.mphi / phi.sphi spec API.
 *
 * FONSE = ROC + nonsense-error selection.  The per-codon log-odds are
 *
 *     eta_i = -dM_i - dEta_i * phi - dOmega_i * phi * beta(pos),
 *     beta(pos) = a1 + a2 * pos,        reference codon: eta = 0,
 *
 * with P(codon i at position pos) = softmax(eta)_i.  dEta is the position-
 * INDEPENDENT elongation selection (identical to ROC's selection term);
 * dOmega is the position-DEPENDENT nonsense-error selection.  As dOmega -> 0
 * FONSE reduces exactly to ROC (roc_sphi_est.stan).
 *
 * Unlike ROC, the per-AA normalising constant Z depends on the codon's
 * position (beta(pos) sits inside the exponent), so the likelihood CANNOT be
 * collapsed to per-gene codon counts.  The data is therefore PER-POSITION
 * (one row per sense codon occurrence): obs_gene / obs_aa / obs_pos /
 * obs_choice.  Rows are grouped by gene (gene_start/gene_end) so reduce_sum
 * can slice over genes.
 *
 * a1 (initiation cost) and a2 (elongation cost) are DATA constants here
 * (default 4 / 4, the native fixed-by-default values).  They form a weakly
 * identified ridge with phi-scale; promote to parameters only with an
 * informative prior.
 *
 * The phi.mphi / phi.sphi spec fields and the centered/non-centered switch are
 * identical to roc_sphi_est.stan -- see that file's header for the full
 * description of every phi_* / sphi_* / noncentered field.
 *
 * THREADING: compile with cpp_options=list(stan_threads=TRUE) and pass
 * threads_per_chain > 1 to enable reduce_sum parallelism.
 * ============================================================================ */

functions {
    real partial_sum(array[] int slice_g, int start, int end,
                          int A,
                          array[] int aa_start, array[] int aa_end,
                          array[] int gene_start, array[] int gene_end,
                          array[] int obs_aa, vector obs_pos,
                          array[] int obs_choice,
                          real a1, real a2,
                          vector dM, vector dEta, vector dOmega,
                          vector phi) {
        real lp = 0;
        int n_slice = size(slice_g);
        for (i in 1:n_slice) {
            int g = slice_g[i];
            real phi_g = phi[g];
            for (r in gene_start[g]:gene_end[g]) {
                int a   = obs_aa[r];
                int s   = aa_start[a];
                int e   = aa_end[a];
                int n   = e - s + 1;                 // number of NON-reference codons
                real beta = a1 + a2 * obs_pos[r];
                // eta_full[1] = reference codon (eta = 0); 2..n+1 = non-reference.
                vector[n + 1] eta_full;
                eta_full[1] = 0;
                for (j in 1:n) {
                    eta_full[j + 1] = -dM[s - 1 + j]
                                      - dEta[s - 1 + j] * phi_g
                                      - dOmega[s - 1 + j] * phi_g * beta;
                }
                lp += eta_full[obs_choice[r]] - log_sum_exp(eta_full);
            }
        }
        return lp;
    }
}

data {
    int<lower=1> G;                       // number of genes
    int<lower=1> A;                       // number of amino acids (>=2-codon)
    int<lower=1> K;                       // number of non-reference codon params
    int<lower=1> N;                       // number of per-position observations
    array[A] int<lower=1> aa_start;       // first param index for each AA
    array[A] int<lower=1> aa_end;         // last  param index for each AA
    array[G] int<lower=1> gene_start;     // first obs row for each gene
    array[G] int<lower=0> gene_end;       // last  obs row for each gene (>= start-1)
    array[N] int<lower=1> obs_aa;         // AA index (1..A) of each observation
    vector[N]             obs_pos;        // codon position (beta = a1 + a2*pos)
    array[N] int<lower=1> obs_choice;     // chosen codon within AA: 1=ref, 2..n+1=nonref

    real a1;                              // initiation cost (data constant; default 4)
    real a2;                              // elongation cost (data constant; default 4)

    vector[K] dM_prior_mean;
    vector<lower=0>[K] dM_prior_sd;
    vector[K] dEta_prior_mean;
    vector<lower=0>[K] dEta_prior_sd;
    vector[K] dOmega_prior_mean;
    vector<lower=0>[K] dOmega_prior_sd;

    // phi.sphi spec
    real<lower=0> sphi_low;
    real          sphi_high;
    int<lower=0, upper=1> sphi_prior_type;
    real          sphi_prior_mean;
    real<lower=0> sphi_prior_sd;

    // phi.mphi spec
    int<lower=0, upper=2> phi_mphi_mode;
    int<lower=0, upper=4> phi_mphi_statistic;
    real<lower=0>         phi_mphi_value;
    real                  phi_mphi_fixed;
    int<lower=0, upper=1> phi_mphi_prior_type;
    real                  phi_mphi_prior_mean;
    real<lower=0>         phi_mphi_prior_sd;

    int<lower=0, upper=1> noncentered;

    int<lower=1> grainsize;
}

transformed data {
    array[G] int gene_indices;
    for (g in 1:G) gene_indices[g] = g;
}

parameters {
    vector[K] dM;
    vector[K] dEta;
    vector[K] dOmega;
    vector[G] latent_phi;     // log_phi (centered) or z_phi (non-centered)
    real<lower=sphi_low, upper=sphi_high> sphi;
    real mphi_param;
}

transformed parameters {
    real s2   = square(sphi);
    real mphi;
    if (phi_mphi_mode == 2) {
        mphi = mphi_param;
    } else if (phi_mphi_mode == 1) {
        mphi = phi_mphi_fixed;
    } else if (phi_mphi_statistic == 0) {
        mphi = log(phi_mphi_value) - 0.5 * s2;
    } else if (phi_mphi_statistic == 1) {
        mphi = log(phi_mphi_value);
    } else if (phi_mphi_statistic == 2) {
        mphi = log(phi_mphi_value) + s2;
    } else if (phi_mphi_statistic == 3) {
        mphi = 0.5 * (log(phi_mphi_value / (exp(s2) - 1.0)) - s2);
    } else {
        mphi = 0.5 * (log(square(phi_mphi_value) / (exp(s2) - 1.0)) - s2);
    }

    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
    vector<lower=0>[G] phi = exp(log_phi);
}

model {
    dM     ~ normal(dM_prior_mean,     dM_prior_sd);
    dEta   ~ normal(dEta_prior_mean,   dEta_prior_sd);
    dOmega ~ normal(dOmega_prior_mean, dOmega_prior_sd);

    if (sphi_prior_type == 1) {
        sphi ~ normal(sphi_prior_mean, sphi_prior_sd);
    }

    if (phi_mphi_mode == 2) {
        if (phi_mphi_prior_type == 1) {
            mphi_param ~ normal(phi_mphi_prior_mean, phi_mphi_prior_sd);
        }
    } else {
        mphi_param ~ std_normal();
    }

    if (noncentered) {
        latent_phi ~ std_normal();
    } else {
        latent_phi ~ normal(mphi, sphi);
    }

    target += reduce_sum(partial_sum, gene_indices, grainsize,
                         A, aa_start, aa_end, gene_start, gene_end,
                         obs_aa, obs_pos, obs_choice,
                         a1, a2, dM, dEta, dOmega, phi);
}
