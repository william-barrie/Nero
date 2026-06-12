#!/usr/bin/env Rscript
# =============================================================================
# f4_common.R
# Shared configuration + data loading + f4 / f4-ratio machinery for the HLA
# African-admixture analysis. Sourced by both:
#   - f4_HLA_analysis_v5.R   (the scan)
#   - f4_sanity_checks.R     (well-known-admixture positive controls)
# so both compute f4 in byte-identical fashion (same streaming frequency
# reader, same block jackknife).
#
# Pure base R (no tidyverse), so it is cheap to source from a diagnostic run.
# See README.md for the memory model and the f4 sign convention.
# =============================================================================

# ---- Run configuration (env-overridable; defaults reproduce v4 behaviour) ---
DATADIR    <- Sys.getenv("NERO_DATADIR", ".")
OUTDIR     <- Sys.getenv("NERO_OUTDIR",  ".")
CHUNK_SNPS <- as.integer(Sys.getenv("NERO_CHUNK_SNPS", "50000"))
if (is.na(CHUNK_SNPS) || CHUNK_SNPS < 1L) CHUNK_SNPS <- 50000L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

out_path <- function(...) file.path(OUTDIR,  ...)
in_path  <- function(...) file.path(DATADIR, ...)

# ---- File prefixes ----------------------------------------------------------
CHR6_PREFIX <- "6.1000g.sea.nocommas"
CHR1_PREFIX <- "1.1000g.sea.nocommas"

# ---- Region & block geometry ------------------------------------------------
HLA_START <- 28000000
HLA_END   <- 33000000
BLOCK_CHR <- 5e6
BLOCK_HLA <- 5e5
BLOCK_WIN <- 5e4
WIN_SIZE  <- 5e5
WIN_STEP  <- 1e5

# =============================================================================
# Data loading: stream .geno -> per-population allele-frequency matrix
# =============================================================================

read_ind <- function(prefix) {
    ind_alt  <- paste0(prefix, ".ind_")
    ind_path <- if (file.exists(ind_alt)) ind_alt else paste0(prefix, ".ind")
    read.table(ind_path, col.names = c("sample", "sex", "pop"),
               stringsAsFactors = FALSE, fill = TRUE)
}

# Read .ind + .snp, then stream .geno in chunks of `chunk_snps` rows, parsing
# ONLY the columns belonging to `want_pops`, and accumulate a per-SNP,
# per-population REF-allele frequency matrix. The full sample-level genotype
# matrix is never materialised. (Verified identical to a dense-matrix read.)
read_eigenstrat_freq <- function(prefix, want_pops, chunk_snps = CHUNK_SNPS) {
    cat(sprintf("  Reading: %s\n", basename(prefix)))

    ind       <- read_ind(prefix)
    n_samples <- nrow(ind)

    snp <- read.table(paste0(prefix, ".snp"),
                      colClasses = c("NULL", "character", "NULL",
                                     "integer", "NULL", "NULL"))
    names(snp) <- c("chr", "pos")
    n_snp <- nrow(snp)

    pop_cols <- lapply(want_pops, function(p) which(ind$pop == p))
    names(pop_cols) <- want_pops
    present  <- want_pops[lengths(pop_cols) > 0]
    if (length(present) == 0L)
        stop("None of the requested populations are present in ", prefix)
    pop_cols <- pop_cols[present]

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
# f4 machinery (frequency-matrix backed)
# =============================================================================

apply_region_mask <- function(dat, start, end) {
    mask <- rep(TRUE, nrow(dat$snp))
    if (!is.null(start)) mask <- mask & (dat$snp$pos >= start)
    if (!is.null(end))   mask <- mask & (dat$snp$pos <= end)
    list(snp = dat$snp[mask, , drop = FALSE],
         freq = dat$freq[mask, , drop = FALSE])
}

have_pops <- function(dat, pops) all(pops %in% colnames(dat$freq))

# f4(A, B; C, D) = mean_SNP (pA - pB)(pC - pD), block (leave-one-out) jackknife
# SE over physical-position blocks of `block_size`.
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
