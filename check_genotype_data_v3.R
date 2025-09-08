# Basic QC reports for each breed (biallelic SNPs only)
# ---------------------------------------------------
# ---- packages ----
need <- c("readr","dplyr","ggplot2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(need, library, character.only = TRUE)

# ---- helpers & thresholds ----
neglog10 <- function(p) -log10(pmax(p, .Machine$double.xmin))  # avoid -Inf

# MAF for biallelic data: use ALT_FREQS
# ALT in plink2 --freq reflects non-reference (VCF/PVAR) by default; it is not guaranteed to be the minor allele.

maf_from_afreq <- function(df) {
  stopifnot("ALT_FREQS" %in% names(df))
  df %>%
    mutate(
      # ensure numeric and clamp to [0,1] just in case
      ALT_FREQS = suppressWarnings(as.numeric(ALT_FREQS)),
      ALT_FREQS = pmin(pmax(ALT_FREQS, 0), 1),
      MAF_calc  = pmin(ALT_FREQS, 1 - ALT_FREQS)
    )
}

# ---- thresholds (tune to taste) ----
mind_thr <- 0.05  # per-sample missingness
geno_thr <- 0.05  # per-variant missingness
hwe_thr  <- 1e-6  # HWE p-value
maf_thr  <- 0.05  # MAF

# ---- breeds & base paths ----
base_dir <- "02_reports"
breeds   <- c("holstein","jersey")

# collect per-breed summary rows here
summary_rows <- list()

# ---- Run for each breed ----
# For holstein
breed = breeds[1]
# For jersey
breed = breeds[2]

# ---- Automatic run across breeds ----
for (breed in breeds) {
  message("=== Processing breed: ", breed, " ===")
  
  # Matches your plink2 --out prefixes:
  #   02_reports/holstein/NIAS_ibv3_296ea.holstein.autosomesX
  #   02_reports/jersey/NIAS_ibv3_296ea.jersey.autosomesX
  prefix <- file.path(base_dir, breed, sprintf("NIAS_ibv3_296ea.%s.autosomesX", breed))
  
  f_smiss <- paste0(prefix, ".smiss")
  f_vmiss <- paste0(prefix, ".vmiss")
  f_hardy <- paste0(prefix, ".hardy")
  f_afreq <- paste0(prefix, ".afreq")
  
  out_dir  <- file.path(base_dir, breed)
  figs_dir <- file.path(out_dir, "figs")
  dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- load ----
  smiss <- readr::read_tsv(f_smiss, show_col_types = FALSE) %>%
    dplyr::rename(FID = `#FID`)  # "#FID" -> FID
  
  vmiss <- readr::read_tsv(f_vmiss, show_col_types = FALSE) %>%
    dplyr::rename(CHROM = `#CHROM`)
  
  hardy <- readr::read_tsv(f_hardy, show_col_types = FALSE) %>%
    dplyr::rename(CHROM = `#CHROM`) %>%
    mutate(minuslog10P = neglog10(P))
  
  afreq <- readr::read_tsv(f_afreq, show_col_types = FALSE) %>%
    dplyr::rename(CHROM = `#CHROM`) %>%
    maf_from_afreq()
  
  # Peek
  message("smiss head:"); print(utils::head(smiss))
  message("vmiss head:"); print(utils::head(vmiss))
  message("hardy head:"); print(utils::head(hardy))
  message("afreq head:"); print(utils::head(afreq))
  
  # ---- plots ----
  p1 <- ggplot(smiss, aes(x = F_MISS)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = mind_thr) +
    labs(title = paste0(breed, ": Sample missingness (F_MISS)"),
         x = "F_MISS", y = "Count")
  
  p2 <- ggplot(vmiss, aes(x = F_MISS)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = geno_thr) +
    labs(title = paste0(breed, ": Variant missingness (F_MISS)"),
         x = "F_MISS", y = "Count")
  
  p3 <- ggplot(hardy, aes(x = minuslog10P)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = -log10(hwe_thr)) +
    labs(title = paste0(breed, ": HWE (-log10 p)"),
         x = "−log10(p)", y = "Count")
  
  p4 <- ggplot(afreq, aes(x = MAF_calc)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = maf_thr) +
    labs(title = paste0(breed, ": MAF distribution"),
         x = "MAF", y = "Count")
  
  print(p1); print(p2); print(p3); print(p4)
  
  ggsave(file.path(figs_dir, "smiss_hist.png"), p1, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figs_dir, "vmiss_hist.png"), p2, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figs_dir, "hardy_hist.png"), p3, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figs_dir, "maf_hist.png"),   p4, width = 6, height = 4, dpi = 300)
  
  # ---- outlier tables (breed-suffixed objects) ----
  assign(paste0("bad_samples_",   breed),
         smiss %>% filter(F_MISS > mind_thr) %>% select(FID, IID, F_MISS),
         envir = .GlobalEnv)
  
  assign(paste0("bad_snps_miss_", breed),
         vmiss %>% filter(F_MISS > geno_thr) %>% select(CHROM, ID, F_MISS),
         envir = .GlobalEnv)
  
  assign(paste0("bad_snps_hwe_",  breed),
         hardy %>% filter(P < hwe_thr) %>% select(CHROM, ID, P),
         envir = .GlobalEnv)
  
  assign(paste0("bad_snps_maf_",  breed),
         afreq %>% filter(MAF_calc < maf_thr) %>% select(CHROM, ID, MAF_calc),
         envir = .GlobalEnv)
  
  # (optional) if you want filenames to include the breed too, even though
  # they already live in per-breed folders:
  readr::write_csv(get(paste0("bad_samples_",   breed)), file.path(out_dir, paste0("bad_samples_over_mind.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_miss_", breed)), file.path(out_dir, paste0("bad_snps_over_geno.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_hwe_",  breed)), file.path(out_dir, paste0("bad_snps_hwe.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_maf_",  breed)), file.path(out_dir, paste0("bad_snps_low_maf.", breed, ".csv")))
  
  
  # ---- combined per-variant table ----
  vmiss <- vmiss %>% mutate(CHROM = as.character(CHROM))
  hardy <- hardy %>% mutate(CHROM = as.character(CHROM))
  afreq <- afreq %>% mutate(CHROM = as.character(CHROM))
  
  variant_qc <- vmiss %>%
    select(CHROM, ID, F_MISS) %>%
    left_join(select(hardy, CHROM, ID, P),        by = c("CHROM","ID")) %>%
    left_join(select(afreq, CHROM, ID, MAF_calc), by = c("CHROM","ID")) %>%
    mutate(
      fail_miss = F_MISS > geno_thr,
      fail_hwe  = P < hwe_thr,
      fail_maf  = MAF_calc < maf_thr,
      fail_any  = fail_miss | fail_hwe | fail_maf
    )
  
  readr::write_csv(variant_qc, file.path(out_dir, paste0("variant_qc_summary.",breed,".csv")))
  
  n_fail_any <- sum(variant_qc$fail_any, na.rm = TRUE)
  
  message(sprintf("[%s] SNPs failing any filter: %d", breed, n_fail_any))
  
  # ---- summary row (compute counts from locals) ----
  samples_over_mind <- sum(smiss$F_MISS > mind_thr, na.rm = TRUE)
  snps_over_geno    <- sum(vmiss$F_MISS > geno_thr, na.rm = TRUE)
  snps_fail_hwe     <- sum(hardy$P < hwe_thr, na.rm = TRUE)
  snps_below_maf    <- sum(afreq$MAF_calc < maf_thr, na.rm = TRUE)
  
  summary_row <- tibble::tibble(
    breed = breed,
    n_samples = nrow(smiss),
    n_variants = nrow(vmiss),
    samples_over_mind = samples_over_mind,
    snps_over_geno = snps_over_geno,
    snps_fail_hwe = snps_fail_hwe,
    snps_below_maf_threshold = snps_below_maf,
    snps_fail_any_filter = n_fail_any
  )
  
  # write per-breed summary and keep in memory
  readr::write_csv(summary_row, file.path(out_dir, "qc_summary_counts.csv"))
  summary_rows[[breed]] <- summary_row
}

# ---- combine summaries across breeds ----
summary_all <- dplyr::bind_rows(summary_rows)
readr::write_csv(summary_all, file.path(base_dir, "qc_summary_by_breed.csv"))

message("Done. Wrote combined summaries to: ",
        normalizePath(file.path(base_dir, "qc_summary_by_breed.csv"), winslash = "/"))



################################################################################
#### Use R markdown!
rmarkdown::render("check_genotype_data.Rmd")
# or with parameters:
rmarkdown::render(
  "qc_biallelic_breed_reports.Rmd",
  params = list(
    base_dir = "02_reports",
    breeds   = c("holstein","jersey"),
    mind_thr = 0.05,
    geno_thr = 0.05,
    hwe_thr  = 1e-6,
    maf_thr  = 0.05
  )
)

rmarkdown::render("check_genotype_data_V2.Rmd")
browseURL("check_genotype_data_V2.html")

rmarkdown::render("check_genotype_data_v3.Rmd")
browseURL("check_genotype_data_V3.html")

pagedown::chrome_print("check_genotype_data_v3.html")

rmarkdown::render("check_genotype_data_v3.Rmd", 
                  output_format = "pagedown::html_paged")

pagedown::chrome_print("check_genotype_data_v3.html", output = "qc_workflow.pdf")

webshot2::rmdshot("check_genotype_data_v3.Rmd", "qc_workflow.pdf")

## Making a book?
install.packages(c("bookdown", "rmarkdown"))
# mkdir my_book
rmarkdown::clean_site()
bookdown::render_book("index.Rmd")
browseURL("_book/index.html")
