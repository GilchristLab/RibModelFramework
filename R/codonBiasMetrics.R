## ============================================================================
## codon_bias_metrics.R
##
## Per-gene codon-bias metrics + null-derivation helpers.
##
## Implements Wright (1990) ENC and Novembre (2002) ENC' for use in:
##   - reference codon selection (ticket #7)
##   - iterative null refinement (ticket #8, future)
##   - audit against SCUO (current default)
##
## Standard genetic code AA family sizes (M, W, X excluded; Serine SPLIT into
## S=TCN and Z=AGY -- see build_aa_codon_map for why):
##   2 codons: C, D, E, F, H, K, N, Q, Y, Z   (10 AAs; Z = AGY Ser)
##   3 codons: I                               (1 AA)
##   4 codons: A, G, P, S, T, V                (6 AAs; S = TCN Ser)
##   6 codons: L, R                            (2 AAs)
## Total: 19 AAs after dropping M/W/X and splitting Ser.
##
## Author: MG + Claude, 2026-05-27.  Generated for issue #7 implementation.
## ============================================================================

## ---- helpers ---------------------------------------------------------------

## Build per-AA codon family map from the standard genetic code.
## Returns a named list: aa -> character vector of synonymous codons (sorted).
##
## Self-contained: does not depend on AnaCoDa or any external package.
## Excludes M (Met), W (Trp) -- single-codon AAs that carry no codon bias info.
## Excludes X (stop codons TAA/TAG/TGA) for the same reason.
##
## Serine is split into TWO families (the AnaCoDa S/Z convention; 19 groups):
##   S = TCN  (TCT/TCC/TCA/TCG, the 4-fold box)
##   Z = AGY  (AGT/AGC,         the 2-fold box)
## These two boxes are mutationally DISCONNECTED -- interconverting TCN<->AGY
## needs TWO substitutions, so no single mutation links them.  Lumping them as
## one 6-fold family (as textbook ENC' does) conflates two independent
## mutation/selection processes; for codon-bias work driving dM/dEta priors the
## split is required.  Leu (TTR/CTN) and Arg (CGN/AGR) are NOT split: their
## boxes ARE linked by single mutations (e.g. TTG<->CTG, CGA<->AGA).
build_aa_codon_map <- function() {
    std_code <- c(
        TTT="F", TTC="F", TTA="L", TTG="L",
        CTT="L", CTC="L", CTA="L", CTG="L",
        ATT="I", ATC="I", ATA="I", ATG="M",
        GTT="V", GTC="V", GTA="V", GTG="V",
        TCT="S", TCC="S", TCA="S", TCG="S",
        CCT="P", CCC="P", CCA="P", CCG="P",
        ACT="T", ACC="T", ACA="T", ACG="T",
        GCT="A", GCC="A", GCA="A", GCG="A",
        TAT="Y", TAC="Y", TAA="X", TAG="X",
        CAT="H", CAC="H", CAA="Q", CAG="Q",
        AAT="N", AAC="N", AAA="K", AAG="K",
        GAT="D", GAC="D", GAA="E", GAG="E",
        TGT="C", TGC="C", TGA="X", TGG="W",
        CGT="R", CGC="R", CGA="R", CGG="R",
        AGT="Z", AGC="Z", AGA="R", AGG="R",   # AGY Ser -> Z (split from TCN Ser)
        GGT="G", GGC="G", GGA="G", GGG="G"
    )
    aas <- sort(setdiff(unique(std_code), c("M", "W", "X")))
    out <- lapply(aas, function(aa) sort(names(std_code)[std_code == aa]))
    names(out) <- aas
    out
}

## Group AAs by codon-family size (returns named list: "2", "3", "4", "6").
group_aa_by_d <- function(aa_codon_map) {
    d <- vapply(aa_codon_map, length, integer(1L))
    split(names(aa_codon_map), as.character(d))
}

## ---- null-derivation -------------------------------------------------------

## Genome-wide per-AA null codon frequencies.
##
## Inputs:
##   codon.counts: G x C matrix of per-gene codon counts (cols named by codon).
##   aa_codon_map: as returned by build_aa_codon_map().
##
## Returns: named list aa -> numeric vector of codon frequencies summing to 1,
## with names = the codons in aa_codon_map[[aa]] (same order).
derive_null_from_genome <- function(codon.counts, aa_codon_map) {
    totals <- colSums(codon.counts)
    out <- lapply(aa_codon_map, function(cs) {
        x <- totals[cs]
        if (sum(x) == 0) return(rep(1/length(cs), length(cs)))   # uniform fallback
        as.numeric(x / sum(x))
    })
    out <- mapply(function(probs, cs) setNames(probs, cs), out, aa_codon_map,
                  SIMPLIFY = FALSE)
    names(out) <- names(aa_codon_map)
    out
}

## Per-AA null codon frequencies from a specified gene subset.
##
## Inputs:
##   codon.counts: G x C matrix.
##   gene_subset:  integer or logical indices into rows of codon.counts.
##   aa_codon_map: as returned by build_aa_codon_map().
##
## Returns: same shape as derive_null_from_genome().
derive_null_from_subset <- function(codon.counts, gene_subset, aa_codon_map) {
    derive_null_from_genome(codon.counts[gene_subset, , drop = FALSE], aa_codon_map)
}

## ---- per-AA F-statistics ---------------------------------------------------

## Wright (1990) unbiased F_a: per-AA homozygosity against UNIFORM null.
##
## F_a = (n * sum(p^2) - 1) / (n - 1)
##
## Floor at 1/d (uniform usage) for stability.  Returns NA if n_a < min_n.
F_wright <- function(counts, min_n = 5L) {
    n <- sum(counts)
    if (n < min_n) return(NA_real_)
    d <- length(counts)
    if (d < 2L) return(NA_real_)
    p <- counts / n
    F_ <- (n * sum(p^2) - 1) / (n - 1)
    if (is.nan(F_)) return(NA_real_)
    max(F_, 1 / d)
}

## Novembre (2002) F'_a: per-AA homozygosity against a GENERAL null.
##
## F'_a = (X^2 + n_a - d) / (d * (n_a - 1))    [Novembre 2002, eq.]
## where X^2 = sum_c (count_c - n_a * p_null_c)^2 / (n_a * p_null_c)
##       d   = number of codons in the AA family
##
## Reduces to Wright's F when p_null is uniform.
## Bounded [1/d, 1].  Returns NA if n_a < min_n.
F_novembre <- function(counts, null_probs, min_n = 5L) {
    n <- sum(counts)
    if (n < min_n) return(NA_real_)
    d <- length(counts)
    if (d < 2L) return(NA_real_)
    if (length(null_probs) != d)
        stop("F_novembre: counts and null_probs have different length")
    expected <- n * null_probs
    ok <- expected > 0
    if (!any(ok)) return(NA_real_)
    X2 <- sum((counts[ok] - expected[ok])^2 / expected[ok])
    F_ <- (X2 + n - d) / (d * (n - 1))
    if (is.nan(F_) || is.infinite(F_)) return(NA_real_)
    min(max(F_, 1 / d), 1)
}

## ---- per-gene ENC and ENC' -------------------------------------------------

## Wright (1990) ENC for a single gene.
##
## ENC = sum_{family size k} N_k / mean(F_a over AAs in that family)
##
## Standard formulation: ENC = 2 + 9/Fbar_2 + 1/Fbar_3 + 5/Fbar_4 + 3/Fbar_6
## (the constant "2" is the contribution from M, W, which we exclude;
##  computed dynamically based on present family sizes).
##
## Returns NA if too many F_a values are NA (defaults: >50% of AAs missing).
.gene_enc_from_F <- function(F_per_aa, aa_codon_map, missing_max_frac = 0.5) {
    n_total <- length(F_per_aa)
    n_missing <- sum(is.na(F_per_aa))
    if (n_missing > missing_max_frac * n_total) return(NA_real_)

    by_d <- group_aa_by_d(aa_codon_map)
    enc <- 0
    for (d_str in names(by_d)) {
        aas_in_family <- by_d[[d_str]]
        F_in_family <- F_per_aa[aas_in_family]
        F_in_family <- F_in_family[!is.na(F_in_family)]
        if (length(F_in_family) == 0L) next
        Fbar <- mean(F_in_family)
        if (Fbar <= 0) Fbar <- 1 / as.integer(d_str)
        enc <- enc + length(aas_in_family) / Fbar
    }
    enc
}

## Per-gene Wright ENC (uniform null), returned as numeric vector of length G.
calc_enc <- function(codon.counts, aa_codon_map, min_n_per_aa = 5L) {
    G <- nrow(codon.counts)
    out <- numeric(G)
    for (g in seq_len(G)) {
        F_g <- vapply(aa_codon_map, function(cs) {
            F_wright(codon.counts[g, cs], min_n = min_n_per_aa)
        }, numeric(1L))
        out[g] <- .gene_enc_from_F(F_g, aa_codon_map)
    }
    names(out) <- rownames(codon.counts)
    out
}

## Per-gene Novembre ENC' (against a supplied per-AA null).
##
## If null_freqs is NULL, uses genome-wide null (derive_null_from_genome on
## the supplied codon.counts).
calc_enc_prime <- function(codon.counts, aa_codon_map,
                           null_freqs = NULL, min_n_per_aa = 5L) {
    if (is.null(null_freqs))
        null_freqs <- derive_null_from_genome(codon.counts, aa_codon_map)
    if (!all(names(aa_codon_map) %in% names(null_freqs)))
        stop("calc_enc_prime: null_freqs missing some AAs")

    G <- nrow(codon.counts)
    out <- numeric(G)
    for (g in seq_len(G)) {
        F_g <- vapply(names(aa_codon_map), function(aa) {
            cs <- aa_codon_map[[aa]]
            F_novembre(codon.counts[g, cs], null_freqs[[aa]], min_n = min_n_per_aa)
        }, numeric(1L))
        out[g] <- .gene_enc_from_F(F_g, aa_codon_map)
    }
    names(out) <- rownames(codon.counts)
    out
}
