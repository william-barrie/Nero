#!/usr/bin/env Rscript
# =============================================================================
# f4_sanity_checks.R
# Positive/negative controls for the f4 set-up: a panel of WELL-KNOWN admixture
# events with textbook expected directions, computed with the SAME machinery as
# the HLA scan (f4_common.R: same streaming reader, same block jackknife).
#
# If the pipeline is working, the positive controls recover the expected sign
# with large |z|, and the null control stays near zero. In particular the
# Africa-into-Southern-Europe control validates exactly the kind of signal the
# HLA analysis is built to detect.
#
# Convention used here (intuitive form):
#     f4(Donor, Outgroup; Test, Control)
#   is POSITIVE when `Test` shares more drift with `Donor` than `Control` does
#   (i.e. Test has more Donor-related ancestry). This is the same compute_f4()
#   used by the scan, just with arguments ordered for readability.
#
# Computed genome-wide on chr1 (the analysis baseline), 5 Mb jackknife blocks.
# Output: sanity_checks.csv + a PASS/FAIL table to stderr.
# =============================================================================

SCRIPT_DIR <- tryCatch({
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
source(file.path(SCRIPT_DIR, "f4_common.R"))

OUTGROUP <- "Siberia_UpperPaleolithic_UstIshim"

# ---- Control panel ----------------------------------------------------------
# expect: "+" positive control (Test > Control in Donor ancestry),
#         "0" null control (no differential ancestry, |z| should stay small).
panel <- rbind.data.frame(
  c("Steppe ancestry higher in N than S Europe (FIN>TSI)",
    "Yamnaya", "FIN", "TSI", "+"),
  c("Steppe ancestry higher in N than S Europe (CEU>IBS)",
    "Yamnaya", "CEU", "IBS", "+"),
  c("CHG ancestry in Yamnaya (Yamnaya = EHG + CHG)",
    "CHG", "Yamnaya", "EHG", "+"),
  c("Hunter-gatherer resurgence in Late vs Early Neolithic",
    "WHG", "FarmerLate", "FarmerEarly", "+"),
  c("Anatolian-farmer ancestry in S Europe vs hunter-gatherers",
    "FarmerAnatolian", "TSI", "WHG", "+"),
  c("African gene flow higher in Iberia than Finland (S>N)",
    "YRI", "IBS", "FIN", "+"),
  c("N-African gene flow higher in Iberia than Finland (Morocco_HG donor)",
    "Morocco_HG", "IBS", "FIN", "+"),
  c("NULL control: CEU vs GBR have ~equal Steppe ancestry",
    "Yamnaya", "CEU", "GBR", "0"),
  stringsAsFactors = FALSE)
names(panel) <- c("test", "donor", "test_pop", "control", "expect")

Z_THRESH <- 3   # |z| considered a clear signal

cat(strrep("=", 78), "\n", sep="")
cat("f4 SANITY CHECKS — well-known admixture events (chr1, genome-wide)\n")
cat("  f4(Donor, Outgroup; Test, Control) > 0  <=>  Test more Donor-like than Control\n")
cat(strrep("=", 78), "\n\n", sep="")

want_pops <- unique(c(OUTGROUP,
                      panel$donor, panel$test_pop, panel$control))

cat("Loading chr1 (frequency-only, controls' populations)...\n")
dat1 <- read_eigenstrat_freq(in_path(CHR1_PREFIX), want_pops)

missing <- setdiff(want_pops, dat1$pops)
if (length(missing))
    cat(sprintf("\nWARNING: populations absent from chr1 (their tests will be NA): %s\n",
                paste(missing, collapse=", ")))

# ---- Run the panel ----------------------------------------------------------
res <- vector("list", nrow(panel))
for (i in seq_len(nrow(panel))) {
    p <- panel[i, ]
    r <- tryCatch(
        compute_f4(dat1, p$donor, OUTGROUP, p$test_pop, p$control,
                   block_size = BLOCK_CHR),
        error = function(e) list(f4=NA, se=NA, z=NA, p=NA, n_snps=0))
    pass <- if (is.na(r$z)) NA
            else if (p$expect == "+") (r$z >  Z_THRESH)
            else                      (abs(r$z) < Z_THRESH)   # null control
    res[[i]] <- data.frame(
        test       = p$test,
        f4_formula = sprintf("f4(%s, O; %s, %s)", p$donor, p$test_pop, p$control),
        expect     = p$expect,
        f4         = r$f4, se = r$se, z = r$z, p = r$p,
        n_snps     = r$n_snps,
        pass       = pass,
        stringsAsFactors = FALSE)
}
res <- do.call(rbind, res)

write.csv(res, out_path("sanity_checks.csv"), row.names = FALSE)

# ---- Report -----------------------------------------------------------------
cat("\n")
cat(sprintf("  %-55s %4s %10s %8s %9s  %s\n",
            "control", "exp", "f4", "z", "p", "result"))
cat("  ", strrep("-", 95), "\n", sep="")
for (i in seq_len(nrow(res))) {
    r <- res[i, ]
    verdict <- if (is.na(r$pass)) "NA (missing pop)"
               else if (r$pass)   "PASS"
               else               "*** CHECK ***"
    cat(sprintf("  %-55s %4s %10.5f %8.2f %9.1e  %s\n",
                substr(r$test, 1, 55), r$expect,
                r$f4, r$z, r$p, verdict))
}
np  <- sum(res$pass, na.rm = TRUE)
nt  <- sum(!is.na(res$pass))
cat(sprintf("\n  %d / %d controls behaved as expected (|z| threshold = %g).\n",
            np, nt, Z_THRESH))
cat("  Saved: ", out_path("sanity_checks.csv"), "\n", sep="")
if (nt > 0 && np < nt)
    cat("  NOTE: a failing control means the f4 set-up is NOT reproducing a\n",
        "        textbook result -- investigate before trusting the HLA scan.\n", sep="")
cat("\nDone.\n")
