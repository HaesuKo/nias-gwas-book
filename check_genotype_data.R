# 품종별 기본 QC 보고서
# ---------------------------------------------------
# ---- 필요한 패키지 설치 및 불러오기----
need <- c("readr","dplyr","ggplot2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(need, library, character.only = TRUE))

# ---- 보조 함수 & 기준값(=필터링 기준) ----
neglog10 <- function(p) -log10(pmax(p, .Machine$double.xmin))  # -Inf 방지

# 이대립성(biallelic) 데이터의 MAF: ALT_FREQS 사용
# plink2 --freq의 ALT는 기본적으로 비참조 대립유전자(.pvar 파일에서 정의된 ALT)이며, 소수 대립유전자(minor allele)가 아닐 수도 있음. 

maf_from_afreq <- function(df) {
  stopifnot("ALT_FREQS" %in% names(df))
  df %>%
    mutate(
      # 숫자형으로 변환 후 [0,1] 범위로 보정
      ALT_FREQS = suppressWarnings(as.numeric(ALT_FREQS)),
      ALT_FREQS = pmin(pmax(ALT_FREQS, 0), 1),
      MAF_calc  = pmin(ALT_FREQS, 1 - ALT_FREQS)
    )
}

# ---- 기준값 (필요 시 조정) ----
mind_thr <- 0.05  # 개체 단위 유전자형 결측률
geno_thr <- 0.05  # 변이 단위 유전자형 결측률
hwe_thr  <- 1e-6  # 하디–바인베르크(HWE) p-값
maf_thr  <- 0.05  # 소수 대립유전자 빈도(MAF)

# ---- 기본 경로 & 품종 구분(서브디렉토리 설정용) ----
base_dir <- "02_reports"
breeds   <- c("holstein","jersey")

# 품종별 요약 행을 저장할 리스트
summary_rows <- list()

# ---- 품종별 실행 ----
# 전체 품종(breeds 벡터: "holstein", "jersey")을 순차적으로 처리한다. 
#   - breeds[1] = "holstein"
#   - breeds[2] = "jersey"
# 아래 for 루프에서 자동으로 두 품종이 반복 실행된다. 

# ---- 모든 품종 자동 실행 ----
for (breed in breeds) {
  message("=== 품종 처리 중: ", breed, " ===")
  
  # plink2 --out 접두어와 일치:
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
  
  # ---- 데이터 불러오기 ----
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
  
  # 미리보기
  message("smiss 데이터 앞부분:"); print(utils::head(smiss))
  message("vmiss 데이터 앞부분:"); print(utils::head(vmiss))
  message("hardy 데이터 앞부분:"); print(utils::head(hardy))
  message("afreq 데이터 앞부분:"); print(utils::head(afreq))
  
  # ---- 그래프 ----
  # ---- 공통 ggplot 테마 ----
  qc_theme <- theme_minimal(base_size = 12) +  # 기본 글씨 크기
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),  # X축 라벨 크게, 두껍게
      axis.title.y = element_text(size = 14, face = "bold"),  # Y축 라벨 크게, 두껍게
      axis.text.x  = element_text(size = 12),                 # 눈금 텍스트 크기
      axis.text.y  = element_text(size = 12),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5) # 가운데 정렬
    )
  
  p1 <- ggplot(smiss, aes(x = F_MISS)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = mind_thr) +
    labs(title = paste0(breed, ": 개체 단위 결측률 (F_MISS)"),
         x = "F_MISS", y = "개수")+qc_theme
  
  p2 <- ggplot(vmiss, aes(x = F_MISS)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = geno_thr) +
    labs(title = paste0(breed, ": 변이 단위 결측률 (F_MISS)"),
         x = "F_MISS", y = "개수")+qc_theme
  
  p3 <- ggplot(hardy, aes(x = minuslog10P)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = -log10(hwe_thr)) +
    labs(title = paste0(breed, ": 하디-바인베르크 검정 (-log10P)"),
         x = "-log10(P)", y = "개수")+qc_theme
  
  p4 <- ggplot(afreq, aes(x = MAF_calc)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = maf_thr) +
    labs(title = paste0(breed, ": MAF 분포"),
         x = "MAF", y = "개수")+qc_theme
  
  print(p1); print(p2); print(p3); print(p4)
  
  ggsave(file.path(figs_dir, "smiss_hist.png"), p1, width = 3, height = 3, dpi = 300)
  ggsave(file.path(figs_dir, "vmiss_hist.png"), p2, width = 3, height = 3, dpi = 300)
  ggsave(file.path(figs_dir, "hardy_hist.png"), p3, width = 3, height = 3, dpi = 300)
  ggsave(file.path(figs_dir, "maf_hist.png"),   p4, width = 3, height = 3, dpi = 300)
  
  # ---- 이상치(outlier) 테이블 (품종명 접미사 포함 객체) ----
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
  
  # (선택 사항) CSV 파일 저장 및 파일 이름에 품종명을 포함시킴(이미 폴더별 구분되어 있음)
  readr::write_csv(get(paste0("bad_samples_",   breed)), file.path(out_dir, paste0("bad_samples_over_mind.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_miss_", breed)), file.path(out_dir, paste0("bad_snps_over_geno.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_hwe_",  breed)), file.path(out_dir, paste0("bad_snps_hwe.", breed, ".csv")))
  readr::write_csv(get(paste0("bad_snps_maf_",  breed)), file.path(out_dir, paste0("bad_snps_low_maf.", breed, ".csv")))
  
  
  # ---- 품종별 변이 QC 리포트 요약 통합 테이블 ----
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
  # 이 유전체 데이터셋에서 기준 미달 SNP 개수
  n_fail_any <- sum(variant_qc$fail_any, na.rm = TRUE)
  # 이 유전체 데이터셋의 총 SNP 수
  n_snps_total <- dplyr::n_distinct(vmiss$ID)
  cat(sprintf("[%s] 총 SNP 수: %s; 기준 실패 SNP 수: %s (%.2f%%)\n\n",
              breed,
              prettyNum(n_snps_total, big.mark = ","),
              prettyNum(n_fail_any,   big.mark = ","),
              100 * n_fail_any / n_snps_total))
  
  # ---- 요약 행 ----
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
    snps_fail_any_filter = n_fail_any,
    pct_fail_any = round(100 * n_fail_any / nrow(vmiss), 3)
  )
  
  # 품종별 QC 리포트 요약 저장 및 메모리에 보관
  readr::write_csv(summary_row, file.path(out_dir, "qc_summary_counts.csv"))
  summary_rows[[breed]] <- summary_row
}

summary_all <- dplyr::bind_rows(summary_rows)
readr::write_csv(summary_all, file.path(base_dir, "qc_summary_by_breed.csv"))