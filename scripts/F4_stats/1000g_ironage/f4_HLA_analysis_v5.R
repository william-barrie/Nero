#!/usr/bin/env Rscript
# =============================================================================
# f4 / f4-ratio analysis: HLA African admixture — v5
# Scan over candidate ancient African source populations
# =============================================================================
# Changes from v4:
#   - MEMORY FIX. v4 read the whole .geno into an n_SNPs x n_samples integer
#     matrix (e.g. 2.7M x 15,391 ~= 169 GB, plus a much larger transient from
#     strsplit() on single characters). v5 instead:
#       * determines the small set of populations the analysis actually needs
#         (outgroup, sister, targets, references, candidate ancient African
#         sources) BEFORE touching the genotypes;
#       * streams the .geno file in bounded chunks, parsing only the columns
#         belonging to those populations (via utf8ToInt, no per-char strings);
#       * collapses straight to a per-population allele-frequency matrix
#         (n_SNPs x n_pops, ~hundreds of MB), which is all the f4 machinery
#         needs. The full sample-level genotype matrix is never materialised.
#     Frequencies are computed identically to v4, so results are unchanged.
#     Peak memory is now governed by NERO_CHUNK_SNPS (default 50,000 SNPs).
#   - REPRODUCIBILITY. Input/output locations are configurable via environment
#     variables (NERO_DATADIR, NERO_OUTDIR, NERO_CHUNK_SNPS) with the previous
#     behaviour as the default. A run manifest (populations_used.csv), the run
#     parameters (run_params.txt) and sessionInfo (sessionInfo.txt) are written
#     so a run can be audited and reproduced.
#   - INTERPRETATION. f4 sign convention is documented inline (see below).
#     Per-pair multiple-testing correction (Benjamini-Hochberg FDR) is added
#     alongside the raw and Bonferroni counts so "how many pairs are really
#     significant" is reported directly.
#
# Inherited from v4:
#   - Outer loop over candidate ancient African source populations (region
#     matches /Africa/i AND groupAge == "Ancient", n_samples >= MIN_SOURCE_N).
#   - Per-source aggregate summary, source ranking, sliding window for the best
#     triple, and the target x reference / source x target / source x reference
#     heatmaps. YRI retained as the modern West African sister for f4-ratio.
#
# f4 SIGN CONVENTION (read this before interpreting outputs)
# -----------------------------------------------------------------------------
#   compute_f4(A, B; C, D) = mean over SNPs of (pA - pB) * (pC - pD),
#   where p is the REF-allele frequency in each population.
#   The contrast tests use  f4(OUTGROUP, source; target, reference).
#   z_diff_chr1 = (f4_HLA - f4_chr1) / se. A NEGATIVE z_diff_chr1 means the HLA
#   region makes the European `target` look MORE like the ancient African
#   `source` (relative to `reference`) than the genome-wide chr1 baseline does
#   -- i.e. localised excess African-source affinity at HLA. Sources are ranked
#   most-negative-first for exactly this reason.
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# ---- Run configuration (env-overridable; defaults reproduce v4 behaviour) ---
DATADIR    <- Sys.getenv("NERO_DATADIR", ".")
OUTDIR     <- Sys.getenv("NERO_OUTDIR",  ".")
CHUNK_SNPS <- as.integer(Sys.getenv("NERO_CHUNK_SNPS", "50000"))
if (is.na(CHUNK_SNPS) || CHUNK_SNPS < 1L) CHUNK_SNPS <- 50000L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

out_path <- function(...) file.path(OUTDIR,  ...)
in_path  <- function(...) file.path(DATADIR, ...)

# ---- File prefixes / metadata path ------------------------------------------
CHR6_PREFIX <- "6.1000g.sea.nocommas"
CHR1_PREFIX <- "1.1000g.sea.nocommas"
META_PATH   <- "sampleInfo.240504_impute_neosea.v07.240627_impute_iadk_cluster.map.ascii.250601_impute_sea_sampleId.tsv"

# ---- Fixed populations ------------------------------------------------------
OUTGROUP <- "Siberia_UpperPaleolithic_UstIshim"
SISTER   <- "YRI"

# Targets are 1000G European populations only — we are testing them as
# recipients of African admixture. We do NOT test ancient European
# populations as targets in this run.
TARGETS <- c("CEU", "GBR", "TSI", "IBS", "FIN")
REFERENCES <- c("WHG", "EHG", "CHG",
                "FarmerAnatolian", "FarmerEarly", "FarmerMiddle", "FarmerLate",
                "Yamnaya")

# ---- Region & block geometry (unchanged from v3) ----------------------------
HLA_START <- 28000000
HLA_END   <- 33000000
BLOCK_CHR <- 5e6
BLOCK_HLA <- 5e5
BLOCK_WIN <- 5e4
WIN_SIZE  <- 5e5
WIN_STEP  <- 1e5

# ---- Source-scan parameters -------------------------------------------------
MIN_SOURCE_N <- 3   # min samples per candidate ancient African source pop

# =============================================================================
# Data loading: stream .geno -> per-population allele-frequency matrix
# =============================================================================

read_ind <- function(prefix) {
    ind_alt  <- paste0(prefix, ".ind_")
    ind_path <- if (file.exists(ind_alt)) ind_alt else paste0(prefix, ".ind")
    read.table(ind_path, col.names = c("sample", "sex", "pop"),
               stringsAsFactors = FALSE, fill = TRUE)
}

# Read the .ind and .snp tables, then stream the .geno file in chunks of
# `chunk_snps` rows, parsing ONLY the columns that belong to `want_pops`, and
# accumulate a per-SNP, per-population REF-allele frequency matrix. The full
# sample-level genotype matrix is never held in memory.
#
# Frequency definition (identical to v4 pop_freq): for a population, the mean
# over its samples of the REF-allele count (0/1/2; missing=9 -> NA), divided
# by 2. SNPs where a population is entirely missing yield NA for that pop.
read_eigenstrat_freq <- function(prefix, want_pops, chunk_snps = CHUNK_SNPS) {
    cat(sprintf("  Reading: %s\n", basename(prefix)))

    ind       <- read_ind(prefix)
    n_samples <- nrow(ind)

    # We only need chromosome + position downstream; skip the rest to save RAM.
    snp <- read.table(paste0(prefix, ".snp"),
                      colClasses = c("NULL", "character", "NULL",
                                     "integer", "NULL", "NULL"))
    names(snp) <- c("chr", "pos")
    n_snp <- nrow(snp)

    # Sample-column indices per requested population (present pops only).
    pop_cols <- lapply(want_pops, function(p) which(ind$pop == p))
    names(pop_cols) <- want_pops
    present  <- want_pops[lengths(pop_cols) > 0]
    if (length(present) == 0L)
        stop("None of the requested populations are present in ", prefix)
    pop_cols <- pop_cols[present]

    # Union of needed sample columns; map each pop's columns into that subset.
    needed         <- sort(unique(unlist(pop_cols, use.names = FALSE)))
    pop_cols_local <- lapply(pop_cols, function(idx) match(idx, needed))

    freq <- matrix(NA_real_, nrow = n_snp, ncol = length(present),
                   dimnames = list(NULL, present))

    con <- file(paste0(prefix, ".geno"), "rt")
    on.exit(close(con))
    row0 <- 0L
    repeat {
        lines <- readLines(con, n = chunk_snps)
        if (length(lines) == 0L) break

        if (row0 == 0L && length(utf8ToInt(lines[1])) != n_samples)
            stop(sprintf(".geno width (%d) != n_samples in .ind (%d) for %s",
                         length(utf8ToInt(lines[1])), n_samples, prefix))

        # Parse only the needed columns. ASCII '0','1','2','9' -> 0,1,2,9.
        m <- t(vapply(lines,
                      function(l) utf8ToInt(l)[needed] - 48L,
                      integer(length(needed)),
                      USE.NAMES = FALSE))
        m[m == 9L] <- NA_integer_

        rng <- (row0 + 1L):(row0 + nrow(m))
        for (p in present) {
            sub <- m[, pop_cols_local[[p]], drop = FALSE]
            freq[rng, p] <- rowMeans(sub, na.rm = TRUE) / 2
        }
        row0 <- row0 + nrow(m)
    }
    if (row0 != n_snp)
        stop(sprintf(".geno rows (%d) != .snp rows (%d) for %s",
                     row0, n_snp, prefix))

    cat(sprintf("    %d SNPs, %d samples, %d pops total | loaded %d pops (%d samples)\n",
                n_snp, n_samples, length(unique(ind$pop)),
                length(present), length(needed)))

    list(ind = ind, snp = snp, freq = freq, pops = present)
}

# =============================================================================
# f4 machinery (frequency-matrix backed; arithmetic unchanged from v4)
# =============================================================================

apply_region_mask <- function(dat, start, end) {
    mask <- rep(TRUE, nrow(dat$snp))
    if (!is.null(start)) mask <- mask & (dat$snp$pos >= start)
    if (!is.null(end))   mask <- mask & (dat$snp$pos <= end)
    list(snp = dat$snp[mask, , drop = FALSE],
         freq = dat$freq[mask, , drop = FALSE])
}

have_pops <- function(dat, pops) all(pops %in% colnames(dat$freq))

compute_f4 <- function(dat, A, B, C, D,
                       start = NULL, end = NULL, block_size = BLOCK_CHR) {
    na_res <- list(f4=NA, se=NA, z=NA, p=NA, n_snps=0, n_blocks=0)
    if (!have_pops(dat, c(A, B, C, D))) return(na_res)
    sub <- apply_region_mask(dat, start, end)
    if (nrow(sub$snp) < 10) return(na_res)

    f4_snp <- (sub$freq[, A] - sub$freq[, B]) * (sub$freq[, C] - sub$freq[, D])
    valid  <- !is.na(f4_snp)
    f4_snp <- f4_snp[valid]
    pos_v  <- sub$snp$pos[valid]
    n <- length(f4_snp)
    if (n < 10) return(list(f4=NA, se=NA, z=NA, p=NA, n_snps=n, n_blocks=0))

    f4_est    <- mean(f4_snp)
    blocks    <- floor(pos_v / block_size)
    block_ids <- unique(blocks)
    n_blocks  <- length(block_ids)
    if (n_blocks < 2)
        return(list(f4=f4_est, se=NA, z=NA, p=NA,
                    n_snps=n, n_blocks=n_blocks))

    f4_jk  <- sapply(block_ids, function(b) mean(f4_snp[blocks != b]))
    pseudo <- n_blocks * f4_est - (n_blocks - 1) * f4_jk
    se <- sqrt(var(pseudo) / n_blocks)
    z  <- f4_est / se
    list(f4=f4_est, se=se, z=z, p=2*pnorm(-abs(z)),
         n_snps=n, n_blocks=n_blocks)
}

# alpha = f4(sister, outgroup; target, reference) /
#         f4(sister, outgroup; source, reference)
compute_f4_ratio <- function(dat, sister, outgroup, target, source, reference,
                             start = NULL, end = NULL, block_size = BLOCK_CHR) {
    na_res <- list(alpha=NA, se_alpha=NA, z_alpha=NA, p_alpha=NA,
                   f4_num=NA, f4_den=NA, n_snps=0, n_blocks=0)
    if (!have_pops(dat, c(sister, outgroup, target, source, reference)))
        return(na_res)
    sub <- apply_region_mask(dat, start, end)
    if (nrow(sub$snp) < 10) return(na_res)

    pA <- sub$freq[, sister]
    pO <- sub$freq[, outgroup]
    pX <- sub$freq[, target]
    pB <- sub$freq[, source]
    pC <- sub$freq[, reference]

    num_snp <- (pA - pO) * (pX - pC)
    den_snp <- (pA - pO) * (pB - pC)

    valid <- !is.na(num_snp) & !is.na(den_snp)
    num_snp <- num_snp[valid]; den_snp <- den_snp[valid]
    pos_v <- sub$snp$pos[valid]
    n <- length(num_snp)
    if (n < 10) { na_res$n_snps <- n; return(na_res) }

    num_est <- mean(num_snp); den_est <- mean(den_snp)
    if (abs(den_est) < 1e-10) return(na_res)
    alpha_est <- num_est / den_est

    blocks    <- floor(pos_v / block_size)
    block_ids <- unique(blocks)
    n_blocks  <- length(block_ids)
    if (n_blocks < 2)
        return(list(alpha=alpha_est, se_alpha=NA, z_alpha=NA, p_alpha=NA,
                    f4_num=num_est, f4_den=den_est,
                    n_snps=n, n_blocks=n_blocks))

    alpha_jk <- sapply(block_ids, function(b) {
        keep  <- blocks != b
        m_num <- mean(num_snp[keep])
        m_den <- mean(den_snp[keep])
        if (abs(m_den) < 1e-10) NA_real_ else m_num / m_den
    })
    alpha_jk <- alpha_jk[!is.na(alpha_jk)]
    n_eff <- length(alpha_jk)
    if (n_eff < 2)
        return(list(alpha=alpha_est, se_alpha=NA, z_alpha=NA, p_alpha=NA,
                    f4_num=num_est, f4_den=den_est,
                    n_snps=n, n_blocks=n_blocks))

    pseudo   <- n_eff * alpha_est - (n_eff - 1) * alpha_jk
    se_alpha <- sqrt(var(pseudo) / n_eff)
    z_alpha  <- alpha_est / se_alpha
    list(alpha=alpha_est, se_alpha=se_alpha, z_alpha=z_alpha,
         p_alpha=2*pnorm(-abs(z_alpha)),
         f4_num=num_est, f4_den=den_est,
         n_snps=n, n_blocks=n_blocks)
}

safe_diff <- function(a, b) {
    if (is.na(a$f4) || is.na(b$f4) || is.na(a$se) || is.na(b$se))
        return(c(diff=NA, se=NA, z=NA, p=NA))
    d  <- a$f4 - b$f4
    sd <- sqrt(a$se^2 + b$se^2)
    z  <- d / sd
    c(diff=d, se=sd, z=z, p=2*pnorm(-abs(z)))
}

test_f4_hla_vs_control <- function(dat6, dat1, A, B, C, D) {
    f4_hla  <- compute_f4(dat6, A, B, C, D,
                          start=HLA_START, end=HLA_END, block_size=BLOCK_HLA)
    f4_rest <- tryCatch({
        snp6 <- dat6$snp
        keep <- (snp6$pos < HLA_START) | (snp6$pos > HLA_END)
        compute_f4(list(snp=snp6[keep, , drop=FALSE],
                        freq=dat6$freq[keep, , drop=FALSE]),
                   A, B, C, D, block_size=BLOCK_CHR)
    }, error = function(e) list(f4=NA, se=NA, n_snps=NA))
    f4_chr1 <- compute_f4(dat1, A, B, C, D, block_size=BLOCK_CHR)

    d_c1 <- safe_diff(f4_hla, f4_chr1)
    d_c6 <- safe_diff(f4_hla, f4_rest)

    data.frame(
        f4_hla=f4_hla$f4, se_hla=f4_hla$se, n_hla=f4_hla$n_snps,
        f4_chr6rest=f4_rest$f4, se_chr6rest=f4_rest$se,
        f4_chr1=f4_chr1$f4, se_chr1=f4_chr1$se, n_chr1=f4_chr1$n_snps,
        diff_chr1=d_c1["diff"], se_diff_chr1=d_c1["se"],
        z_diff_chr1=d_c1["z"],  p_diff_chr1=d_c1["p"],
        diff_chr6rest=d_c6["diff"], se_diff_chr6rest=d_c6["se"],
        z_diff_chr6rest=d_c6["z"],  p_diff_chr6rest=d_c6["p"],
        row.names=NULL, stringsAsFactors=FALSE)
}

# =============================================================================
# Candidate ancient African source populations
# =============================================================================

get_african_sources <- function(meta_path, ind, min_n) {
    meta <- read.table(meta_path, sep="\t", header=TRUE, quote="",
                       comment.char="", stringsAsFactors=FALSE,
                       na.strings=c("", "NA"))
    afr_pops <- meta %>%
        filter(grepl("Africa", region, ignore.case=TRUE),
               groupAge == "Ancient",
               !is.na(groupLabel)) %>%
        pull(groupLabel) %>% unique()

    counts <- table(ind$pop)
    keep   <- intersect(afr_pops, names(counts))
    keep   <- keep[counts[keep] >= min_n]
    sort(keep)
}

# =============================================================================
# MAIN
# =============================================================================

cat(strrep("=", 70), "\n", sep="")
cat("f4 ANALYSIS: BEST ANCIENT AFRICAN SOURCE FOR HLA SIGNAL — v5\n")
cat(strrep("=", 70), "\n\n", sep="")

cat(sprintf("DATADIR=%s | OUTDIR=%s | CHUNK_SNPS=%d\n\n",
            DATADIR, OUTDIR, CHUNK_SNPS))

# Resolve which populations we need BEFORE loading any genotypes. The candidate
# African sources are derived from the chr6 .ind labels + metadata, which are
# small, so we read those first.
chr6_pre <- in_path(CHR6_PREFIX)
chr1_pre <- in_path(CHR1_PREFIX)

cat("Resolving populations to load...\n")
ind6  <- read_ind(chr6_pre)
pops6 <- unique(ind6$pop)

missing_essential <- setdiff(c(OUTGROUP, SISTER), pops6)
if (length(missing_essential) > 0)
    stop(sprintf("Required populations missing from chr6 data: %s",
                 paste(missing_essential, collapse=", ")))

targets_ok <- TARGETS[TARGETS %in% pops6]
refs_ok    <- REFERENCES[REFERENCES %in% pops6]
sources_ok <- get_african_sources(in_path(META_PATH), ind6, MIN_SOURCE_N)

want_pops <- unique(c(OUTGROUP, SISTER, targets_ok, refs_ok, sources_ok))
cat(sprintf("Populations required by the analysis: %d\n\n", length(want_pops)))

cat("Loading data (streaming, frequency-only)...\n")
dat6 <- read_eigenstrat_freq(chr6_pre, want_pops)
dat1 <- read_eigenstrat_freq(chr1_pre, want_pops)

cat(sprintf("\nTargets present:    %d / %d\n", length(targets_ok), length(TARGETS)))
cat(sprintf("References present: %d / %d\n", length(refs_ok),    length(REFERENCES)))
miss_t <- setdiff(TARGETS,    pops6); if (length(miss_t)) cat("  Missing targets:   ", paste(miss_t, collapse=", "), "\n")
miss_r <- setdiff(REFERENCES, pops6); if (length(miss_r)) cat("  Missing references:", paste(miss_r, collapse=", "), "\n")

cat(sprintf("\nCandidate African sources (n>=%d, present in chr6 .ind): %d\n",
            MIN_SOURCE_N, length(sources_ok)))
counts6 <- table(ind6$pop)
for (s in sources_ok)
    cat(sprintf("    %-45s n=%d\n", s, counts6[s]))

if (length(sources_ok) == 0) stop("No candidate African sources found.")

# Sister/source overlap check (warn only)
sister_ids <- ind6$sample[ind6$pop == SISTER]
for (src in sources_ok) {
    src_ids <- ind6$sample[ind6$pop == src]
    ov <- intersect(sister_ids, src_ids)
    if (length(ov) > 0)
        warning(sprintf("%d sample(s) shared between SISTER (%s) and source (%s)",
                        length(ov), SISTER, src))
}

# ---- Reproducibility: write a manifest of populations actually loaded -------
role_of <- function(p) {
    if (p == OUTGROUP) "outgroup"
    else if (p == SISTER) "sister"
    else if (p %in% targets_ok) "target"
    else if (p %in% refs_ok) "reference"
    else if (p %in% sources_ok) "source"
    else "other"
}
n1 <- table(dat1$ind$pop)
manifest <- tibble(
    population = want_pops,
    role       = vapply(want_pops, role_of, character(1)),
    n_chr6     = as.integer(counts6[want_pops]),
    n_chr1     = as.integer(n1[want_pops]),
    in_chr6    = want_pops %in% dat6$pops,
    in_chr1    = want_pops %in% dat1$pops) %>%
    arrange(role, population)
write.csv(manifest, out_path("populations_used.csv"), row.names=FALSE)
cat(sprintf("Saved: %s\n", out_path("populations_used.csv")))

n_pairs_per_src <- length(refs_ok) * length(targets_ok)
cat(sprintf("\nPairs per source: %d   |   total f4 calls ~ %d\n",
            n_pairs_per_src, n_pairs_per_src * length(sources_ok) * 3))

# -----------------------------------------------------------------------------
# TEST 1 + TEST 2: f4 contrast and f4-ratio across all (source, target, ref)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep="")
cat("TEST 1+2: f4 CONTRAST and f4-RATIO (all source × target × reference)\n")
cat(strrep("=", 70), "\n", sep="")

f4_results    <- list()
ratio_results <- list()
t0 <- Sys.time()

for (si in seq_along(sources_ok)) {
    src <- sources_ok[si]
    cat(sprintf("\n[%s] [%2d/%2d] Source: %s\n",
                format(Sys.time(), "%H:%M:%S"), si, length(sources_ok), src))

    f4_block    <- vector("list", n_pairs_per_src); fk <- 0
    ratio_block <- vector("list", n_pairs_per_src); rk <- 0

    for (ref in refs_ok) {
        for (target in targets_ok) {

            # f4 contrast
            res <- tryCatch(
                test_f4_hla_vs_control(dat6, dat1, OUTGROUP, src, target, ref),
                error = function(e) NULL)
            if (!is.null(res)) {
                res$source <- src; res$target <- target; res$reference <- ref
                fk <- fk + 1; f4_block[[fk]] <- res
            }

            # f4-ratio (HLA vs chr1 baseline)
            r_hla <- tryCatch(
                compute_f4_ratio(dat6, SISTER, OUTGROUP, target, src, ref,
                                 start=HLA_START, end=HLA_END, block_size=BLOCK_HLA),
                error = function(e) NULL)
            r_chr1 <- tryCatch(
                compute_f4_ratio(dat1, SISTER, OUTGROUP, target, src, ref,
                                 block_size=BLOCK_CHR),
                error = function(e) NULL)
            if (!is.null(r_hla) && !is.null(r_chr1) &&
                !is.na(r_hla$alpha) && !is.na(r_chr1$alpha) &&
                !is.na(r_hla$se_alpha) && !is.na(r_chr1$se_alpha)) {
                d  <- r_hla$alpha - r_chr1$alpha
                sd <- sqrt(r_hla$se_alpha^2 + r_chr1$se_alpha^2)
                z  <- d / sd
                rk <- rk + 1
                ratio_block[[rk]] <- data.frame(
                    source=src, target=target, reference=ref,
                    alpha_hla=r_hla$alpha,   se_alpha_hla=r_hla$se_alpha,
                    alpha_chr1=r_chr1$alpha, se_alpha_chr1=r_chr1$se_alpha,
                    delta_alpha=d, se_delta=sd, z_delta=z, p_delta=2*pnorm(-abs(z)),
                    stringsAsFactors=FALSE)
            }
        }
    }

    f4_results[[src]]    <- if (fk > 0) do.call(rbind, f4_block[seq_len(fk)])    else NULL
    ratio_results[[src]] <- if (rk > 0) do.call(rbind, ratio_block[seq_len(rk)]) else NULL

    fb <- f4_results[[src]]
    if (!is.null(fb))
        cat(sprintf("  -> %d pairs | mean z_diff_chr1 = %.2f | min p_diff_chr1 = %.2g\n",
                    nrow(fb),
                    mean(fb$z_diff_chr1, na.rm=TRUE),
                    min(fb$p_diff_chr1,  na.rm=TRUE)))
}

f4_results    <- bind_rows(f4_results)
ratio_results <- bind_rows(ratio_results)

# Multiple-testing correction across all (source, target, reference) pairs.
if (nrow(f4_results) > 0)
    f4_results$fdr_diff_chr1 <- p.adjust(f4_results$p_diff_chr1, method="BH")
if (nrow(ratio_results) > 0)
    ratio_results$fdr_delta  <- p.adjust(ratio_results$p_delta, method="BH")

cat(sprintf("\nElapsed: %s\n", format(Sys.time() - t0)))

write.csv(f4_results,    out_path("f4_contrast_results.csv"), row.names=FALSE)
write.csv(ratio_results, out_path("f4_ratio_results.csv"),    row.names=FALSE)
cat("Saved: f4_contrast_results.csv, f4_ratio_results.csv\n")

# -----------------------------------------------------------------------------
# Per-source summary
# -----------------------------------------------------------------------------
n_total <- nrow(f4_results)

f4_summary <- f4_results %>%
    group_by(source) %>%
    summarise(
        n_pairs              = n(),
        mean_z_diff_chr1     = mean(z_diff_chr1,     na.rm=TRUE),
        mean_z_diff_chr6rest = mean(z_diff_chr6rest, na.rm=TRUE),
        min_p_diff_chr1      = suppressWarnings(min(p_diff_chr1, na.rm=TRUE)),
        n_sig_chr1_05        = sum(p_diff_chr1     < 0.05,           na.rm=TRUE),
        n_sig_chr1_bonf      = sum(p_diff_chr1     < 0.05 / n_total, na.rm=TRUE),
        n_sig_chr1_fdr       = sum(fdr_diff_chr1   < 0.05,           na.rm=TRUE),
        n_sig_chr6rest_05    = sum(p_diff_chr6rest < 0.05,           na.rm=TRUE),
        .groups = "drop")

ratio_summary <- ratio_results %>%
    group_by(source) %>%
    summarise(
        mean_alpha_hla   = mean(alpha_hla,   na.rm=TRUE),
        mean_alpha_chr1  = mean(alpha_chr1,  na.rm=TRUE),
        mean_delta_alpha = mean(delta_alpha, na.rm=TRUE),
        mean_z_delta     = mean(z_delta,     na.rm=TRUE),
        n_sig_delta_05   = sum(p_delta < 0.05, na.rm=TRUE),
        n_sig_delta_fdr  = sum(fdr_delta < 0.05, na.rm=TRUE),
        .groups = "drop")

src_summary <- f4_summary %>%
    left_join(ratio_summary, by="source") %>%
    arrange(mean_z_diff_chr1)
write.csv(src_summary, out_path("source_summary.csv"), row.names=FALSE)

cat("\n--- Source ranking by mean z_diff_chr1 (most negative = strongest signal) ---\n")
cat(sprintf("  %-45s %5s %9s %9s %12s %10s %12s\n",
            "source", "n", "mean_z1", "mean_z6r", "min_p_chr1", "n_sig05", "mean_d_alpha"))
for (i in seq_len(nrow(src_summary))) {
    r <- src_summary[i, ]
    cat(sprintf("  %-45s %5d %9.2f %9.2f %12.2e %10d %12.4f\n",
                r$source, r$n_pairs,
                r$mean_z_diff_chr1, r$mean_z_diff_chr6rest,
                r$min_p_diff_chr1, r$n_sig_chr1_05,
                ifelse(is.na(r$mean_delta_alpha), NA_real_, r$mean_delta_alpha)))
}

# -----------------------------------------------------------------------------
# TEST 3: sliding window for the BEST (source, target, reference) triple
# -----------------------------------------------------------------------------
best_source <- src_summary$source[1]
best_block  <- f4_results %>%
    filter(source == best_source) %>%
    arrange(z_diff_chr1)
best_target <- best_block$target[1]
best_ref    <- best_block$reference[1]

cat("\n", strrep("=", 70), "\n", sep="")
cat("TEST 3: SLIDING WINDOW\n")
cat(strrep("=", 70), "\n", sep="")
cat(sprintf("Best triple: source=%s | target=%s | reference=%s | z_chr1=%.2f\n\n",
            best_source, best_target, best_ref, best_block$z_diff_chr1[1]))

starts <- seq(min(dat6$snp$pos),
              max(dat6$snp$pos) - WIN_SIZE, by=WIN_STEP)

# f4 sliding window
win_f4 <- map_dfr(starts, function(s) {
    e <- s + WIN_SIZE
    res <- tryCatch(
        compute_f4(dat6, OUTGROUP, best_source, best_target, best_ref,
                   start=s, end=e, block_size=BLOCK_WIN),
        error = function(err) NULL)
    if (!is.null(res) && !is.na(res$f4) && res$n_snps >= 10)
        tibble(midpoint=(s+e)/2, f4=res$f4, se=res$se, z=res$z, n_snps=res$n_snps)
    else NULL
})

if (nrow(win_f4) > 0) {
    write.csv(win_f4, out_path("chr6_f4_sliding_window.csv"), row.names=FALSE)
    f4_baseline <- compute_f4(dat1, OUTGROUP, best_source, best_target, best_ref,
                              block_size=BLOCK_CHR)
    y_rng   <- range(win_f4$f4, na.rm=TRUE)
    label_y <- y_rng[1] - 0.05 * diff(y_rng)

    p1 <- ggplot(win_f4, aes(x=midpoint/1e6, y=f4)) +
        annotate("rect", xmin=HLA_START/1e6, xmax=HLA_END/1e6,
                 ymin=-Inf, ymax=Inf, fill="blue", alpha=0.08) +
        geom_ribbon(aes(ymin=f4-1.96*se, ymax=f4+1.96*se),
                    alpha=0.18, fill="steelblue") +
        geom_line(color="grey25", linewidth=0.4) +
        geom_hline(yintercept=0, linetype="dotted", color="grey60") +
        {if (!is.na(f4_baseline$f4))
            geom_hline(yintercept=f4_baseline$f4,
                       linetype="dashed", color="red", linewidth=0.5)} +
        annotate("text", x=(HLA_START+HLA_END)/2e6, y=label_y,
                 label="HLA", color="blue", size=4, fontface="bold") +
        labs(x="Chr6 position (Mb)",
             y=sprintf("f4(%s, %s; %s, %s)",
                       OUTGROUP, best_source, best_target, best_ref),
             title=sprintf("Sliding-window f4 — best source: %s", best_source),
             subtitle=sprintf("target=%s, reference=%s | red dashed = chr1 baseline | %g kb windows, %g kb step",
                              best_target, best_ref, WIN_SIZE/1e3, WIN_STEP/1e3)) +
        theme_minimal(base_size=11) +
        theme(plot.title=element_text(face="bold"))
    ggsave(out_path("chr6_f4_sliding_window.pdf"), p1, width=12, height=4.5)
    cat("Saved: chr6_f4_sliding_window.pdf\n")

    # alpha sliding window
    win_alpha <- map_dfr(starts, function(s) {
        e <- s + WIN_SIZE
        r <- tryCatch(
            compute_f4_ratio(dat6, SISTER, OUTGROUP, best_target, best_source, best_ref,
                             start=s, end=e, block_size=BLOCK_WIN),
            error = function(err) NULL)
        if (!is.null(r) && !is.na(r$alpha) && !is.na(r$se_alpha) && r$n_snps >= 10)
            tibble(midpoint=(s+e)/2, alpha=r$alpha, se=r$se_alpha, n_snps=r$n_snps)
        else NULL
    })

    if (nrow(win_alpha) > 0) {
        write.csv(win_alpha, out_path("chr6_alpha_sliding_window.csv"), row.names=FALSE)
        alpha_baseline <- compute_f4_ratio(dat1, SISTER, OUTGROUP,
                                           best_target, best_source, best_ref,
                                           block_size=BLOCK_CHR)
        a_rng  <- range(win_alpha$alpha, na.rm=TRUE)
        lab_ya <- a_rng[2] + 0.05 * diff(a_rng)
        p2 <- ggplot(win_alpha, aes(x=midpoint/1e6, y=alpha)) +
            annotate("rect", xmin=HLA_START/1e6, xmax=HLA_END/1e6,
                     ymin=-Inf, ymax=Inf, fill="blue", alpha=0.08) +
            geom_ribbon(aes(ymin=alpha-1.96*se, ymax=alpha+1.96*se),
                        alpha=0.18, fill="darkgreen") +
            geom_line(color="grey25", linewidth=0.4) +
            geom_hline(yintercept=0, linetype="dotted", color="grey60") +
            {if (!is.na(alpha_baseline$alpha))
                geom_hline(yintercept=alpha_baseline$alpha,
                           linetype="dashed", color="red", linewidth=0.5)} +
            annotate("text", x=(HLA_START+HLA_END)/2e6, y=lab_ya,
                     label="HLA", color="blue", size=4, fontface="bold") +
            labs(x="Chr6 position (Mb)",
                 y=expression(alpha~"(admixture proportion via best source)"),
                 title=sprintf("f4-ratio alpha — best source: %s", best_source),
                 subtitle=sprintf("target=%s, reference=%s | red dashed = chr1 baseline",
                                  best_target, best_ref)) +
            theme_minimal(base_size=11) +
            theme(plot.title=element_text(face="bold"))
        ggsave(out_path("chr6_alpha_sliding_window.pdf"), p2, width=12, height=4.5)
        cat("Saved: chr6_alpha_sliding_window.pdf\n")
    }
}

# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------

# (a) Original-style target × reference heatmap, for the best source only
best_df <- f4_results %>%
    filter(source == best_source) %>%
    select(target, reference, z_diff_chr1) %>%
    filter(!is.na(target), !is.na(reference))
if (nrow(best_df) > 1) {
    target_order <- best_df %>%
        group_by(target) %>%
        summarise(m = mean(z_diff_chr1, na.rm=TRUE), .groups="drop") %>%
        arrange(m) %>% pull(target)
    best_df$target <- factor(best_df$target, levels=target_order)

    p_best <- ggplot(best_df, aes(x=reference, y=target, fill=z_diff_chr1)) +
        geom_tile(color="white", linewidth=0.25) +
        scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",
                             midpoint=0, na.value="grey85",
                             name=expression(z[diff])) +
        labs(title=sprintf("f4 HLA excess over chr1 baseline — source: %s", best_source),
             subtitle="Blue = HLA more source-like than chr1",
             x="Reference", y="Target") +
        theme_minimal(base_size=10) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              panel.grid=element_blank(),
              plot.title=element_text(face="bold"))
    ggsave(out_path("f4_hla_excess_heatmap_bestsource.pdf"), p_best,
           width=8, height=max(4, length(unique(best_df$target)) * 0.4))
    cat("Saved: f4_hla_excess_heatmap_bestsource.pdf\n")
}

# (b) source × target summary (mean z over references)
src_tgt <- f4_results %>%
    group_by(source, target) %>%
    summarise(mean_z = mean(z_diff_chr1, na.rm=TRUE), .groups="drop")
if (nrow(src_tgt) > 1) {
    src_order <- src_summary %>% pull(source)   # already sorted by mean_z_diff_chr1
    src_tgt$source <- factor(src_tgt$source, levels=src_order)

    p_st <- ggplot(src_tgt, aes(x=target, y=source, fill=mean_z)) +
        geom_tile(color="white", linewidth=0.25) +
        scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",
                             midpoint=0, na.value="grey85",
                             name="mean z\n(over refs)") +
        labs(title="Source × target heatmap (mean z_diff_chr1 over references)",
             subtitle="Blue = HLA more source-like than chr1",
             x="Target", y="Source (sorted by overall signal)") +
        theme_minimal(base_size=9) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              panel.grid=element_blank(),
              plot.title=element_text(face="bold"))
    ggsave(out_path("source_x_target_heatmap.pdf"), p_st,
           width=8, height=max(5, length(src_order) * 0.28),
           limitsize=FALSE)
    cat("Saved: source_x_target_heatmap.pdf\n")
}

# (c) source × reference summary (mean z over targets)
src_ref <- f4_results %>%
    group_by(source, reference) %>%
    summarise(mean_z = mean(z_diff_chr1, na.rm=TRUE), .groups="drop")
if (nrow(src_ref) > 1) {
    src_order <- src_summary %>% pull(source)
    src_ref$source <- factor(src_ref$source, levels=src_order)

    p_sr <- ggplot(src_ref, aes(x=reference, y=source, fill=mean_z)) +
        geom_tile(color="white", linewidth=0.25) +
        scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",
                             midpoint=0, na.value="grey85",
                             name="mean z\n(over targets)") +
        labs(title="Source × reference heatmap (mean z_diff_chr1 over targets)",
             x="Reference", y="Source (sorted by overall signal)") +
        theme_minimal(base_size=9) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              panel.grid=element_blank(),
              plot.title=element_text(face="bold"))
    ggsave(out_path("source_x_reference_heatmap.pdf"), p_sr,
           width=8, height=max(5, length(src_order) * 0.28),
           limitsize=FALSE)
    cat("Saved: source_x_reference_heatmap.pdf\n")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep="")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n", sep="")
cat(sprintf("Sources tested:           %d\n", length(sources_ok)))
cat(sprintf("Pairs per source:         %d\n", n_pairs_per_src))
cat(sprintf("Total pair tests:         %d\n", n_total))
cat(sprintf("Best source (mean z):     %s  (mean_z_diff_chr1 = %.2f)\n",
            best_source, src_summary$mean_z_diff_chr1[1]))
cat(sprintf("Best triple sliding-win:  %s × %s × %s\n",
            best_source, best_target, best_ref))

# ---- Reproducibility: record run parameters and the R session ---------------
writeLines(c(
    "f4_HLA_analysis_v5.R run parameters",
    paste0("date            : ", format(Sys.time())),
    paste0("DATADIR         : ", normalizePath(DATADIR, mustWork=FALSE)),
    paste0("OUTDIR          : ", normalizePath(OUTDIR,  mustWork=FALSE)),
    paste0("CHUNK_SNPS      : ", CHUNK_SNPS),
    paste0("OUTGROUP        : ", OUTGROUP),
    paste0("SISTER          : ", SISTER),
    paste0("TARGETS         : ", paste(TARGETS,    collapse=", ")),
    paste0("REFERENCES      : ", paste(REFERENCES, collapse=", ")),
    paste0("MIN_SOURCE_N    : ", MIN_SOURCE_N),
    paste0("HLA window      : ", HLA_START, "-", HLA_END),
    paste0("sources_ok      : ", paste(sources_ok, collapse=", ")),
    paste0("best_triple     : ", best_source, " / ", best_target, " / ", best_ref)
), out_path("run_params.txt"))
writeLines(capture.output(sessionInfo()), out_path("sessionInfo.txt"))

cat("\nOutputs (in ", normalizePath(OUTDIR, mustWork=FALSE), "):\n", sep="")
cat("  populations_used.csv\n")
cat("  f4_contrast_results.csv\n")
cat("  f4_ratio_results.csv\n")
cat("  source_summary.csv\n")
cat("  chr6_f4_sliding_window.{csv,pdf}\n")
cat("  chr6_alpha_sliding_window.{csv,pdf}\n")
cat("  f4_hla_excess_heatmap_bestsource.pdf\n")
cat("  source_x_target_heatmap.pdf\n")
cat("  source_x_reference_heatmap.pdf\n")
cat("  run_params.txt, sessionInfo.txt\n")
cat("\nDone.\n")
