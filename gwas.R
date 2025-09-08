## ====================================================================
## 패키지 설치 & 로딩 (한 번만 실행하면 됨)
## ====================================================================
cran_pkgs <- c("genio","qqman","sommer","lme4breeding","dplyr","CMplot","stringr")
bioc_pkgs <- c("GENESIS","GWASTools","SNPRelate","BiocParallel","rtracklayer","GenomicRanges")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install_if_missing <- function(pkgs, installer) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) installer(need, ask = FALSE)
}

install_if_missing(cran_pkgs, install.packages)
install_if_missing(bioc_pkgs, BiocManager::install)

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# 단계별 튜토리얼: GWAS를 R에서 수행하기
# 이 튜토리얼은 초보자가 R에서 PLINK 유전자형 데이터 파일과 표현형 데이터를 이용해
# 전장유전체연관분석(GWAS)을 수행하는 과정을 안내한다.
# 1) 필요한 패키지와 데이터를 불러온다.
# 2) SNP chip 유전자형/표현형 데이터를 준비한다.
# 3) 유전체 관계 행렬(GRM)을 계산한다.
# 4) GRM을 사용해 유전력(heritability)을 추정한다.
# 5) GWASTools를 이용해 GWAS용 데이터 구조를 준비한다.
# 6) 개체 간 유전적 관계/집단 구조(GRM 기반 임의효과), 고정효과, 공변량을 반영한 귀무 모형(null model)을 적합한다.
# 7) 유전자형과 표현형 간 연관성 검정을 수행한다.
# 8) p-값을 보정하고 유의한 SNP를 선별한다.
# 9) 유의 SNP에 대해 유전자를 주석(annotation)한다.

# ---- 라이브러리 로드 -------------------------------------------------
library(genio)         # PLINK I/O                [genio]
library(GENESIS)       # 귀무모형/연관성검정        [GENESIS]
library(GWASTools)     # 유전자형 데이터구조/iterator [GWASTools]
library(qqman)         # QQ/맨해튼 플롯           [qqman]
library(lme4breeding)  # LMM, 임의효과 설명력     [lme4breeding]
library(sommer)        # A.mat(GRM)               [sommer]
library(dplyr)         # 데이터 가공              [dplyr]

# ---- Step 2: 유전자형 데이터(PLINK) 불러오기 --------------------------------

bw8_input_corrected <- read_plink('bw8_corrected') # Load PLINK data; assumes files are in working directory.
X_bw8_corrected <- bw8_input_corrected$X # Extract genotype matrix (rows: SNPs, columns: samples, values: 0/1/2).
bim_bw8_corrected <- bw8_input_corrected$bim # Extract BIM data frame (SNP details: chr, id, pos, alt, ref).
fam_bw8_corrected <- bw8_input_corrected$fam # Extract FAM data frame (sample details).

# ---- Step 3: 표현형 데이터 준비 ---------------------------------------
sample_ids   <- colnames(X_bw8_corrected)
mydat        <- merged5[merged5$F2_id %in% sample_ids, ]
utils::head(mydat, 6)                                                     # [utils::head]
mydat$scanID <- rownames(mydat)
utils::str(mydat)                                                         # [utils::str]

# ---- Step 4: GRM 계산 --------------------------------------------------
common_samples  <- base::intersect(colnames(X_bw8_corrected), mydat$F2_id)
X_bw8_subset    <- X_bw8_corrected[, common_samples]
mydat_ordered   <- mydat[base::match(common_samples, mydat$F2_id), ]
base::all.equal(mydat_ordered$F2_id, colnames(X_bw8_subset))              # [base::all.equal]

geno        <- base::t(X_bw8_subset)
geno_coded  <- geno - 1
utils::str(geno_coded)

myKI <- sommer::A.mat(geno_coded)                                         # [sommer::A.mat]
myKI <- myKI + base::diag(1e-4, ncol(myKI), ncol(myKI))                   # [base::diag]

myKI[1:5, 1:5]
colfunc <- grDevices::colorRampPalette(c("steelblue4","springgreen","yellow")) # [grDevices::colorRampPalette]
grDevices::png("bw8_Amatrix_updated.png", width=10, height=6, units="in", res=300) # [grDevices::png]
stats::heatmap(myKI, col = colfunc(100), Colv = "Rowv", symmetric = TRUE) # [stats::heatmap]
grDevices::dev.off()                                                      # [grDevices::dev.off]

# ---- Step 5: GWASTools 데이터 구조 준비 --------------------------------
bim_bw8_corrected <- base::cbind(index = 1:nrow(bim_bw8_corrected), bim_bw8_corrected)
base::all.equal(rownames(X_bw8_subset), bim_bw8_corrected$id)
rownames(X_bw8_subset) <- bim_bw8_corrected$index

dplyr::distinct(bim_bw8_corrected, chr) %>% print()

geno <- GWASTools::MatrixGenotypeReader(                                  # [GWASTools::MatrixGenotypeReader]
  genotype   = X_bw8_subset,
  snpID      = bim_bw8_corrected$index,
  chromosome = base::as.integer(bim_bw8_corrected$chr),
  position   = base::as.integer(bim_bw8_corrected$pos),
  scanID     = base::as.integer(mydat_ordered$scanID),
  autosomeCode = 1L:33L,
  XchromCode   = 34L,
  YchromCode   = 35L,
  XYchromCode  = 36L,
  MchromCode   = 37L
)
genoData <- GWASTools::GenotypeData(geno)                                 # [GWASTools::GenotypeData]
scanAnno <- GWASTools::ScanAnnotationDataFrame(mydat_ordered)             # [GWASTools::ScanAnnotationDataFrame]
scanAnno

# ---- Step 6: 변수 정렬 확인 & LMM 적합 ----------------------------------
rownames(myKI); colnames(myKI)
base::all.equal(rownames(myKI), colnames(myKI))
base::all.equal(colnames(myKI), colnames(X_bw8_subset))
base::all.equal(mydat_ordered$F2_id, colnames(myKI))
rownames(myKI) <- mydat_ordered$scanID
colnames(myKI) <- mydat_ordered$scanID

mix1 <- lme4breeding::lmebreed(                                           # [lme4breeding::lmebreed]
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_sire),
  relmat  = list(scanID = myKI),
  verbose = TRUE,
  data    = mydat_ordered
)
utils::summary(mix1)                                                      # [utils::summary]
vc <- lme4breeding::VarCorr(mix1)                                         # [lme4breeding::VarCorr]
print(vc, comp = c("Variance"))
ve <- attr(lme4breeding::VarCorr(mix1), "sc")^2
h2 <- vc$scanID / (vc$scanID + vc$F1_sire + ve)
base::as.numeric(h2)
vc$F1_sire / (vc$scanID + vc$F1_sire + ve)

mix2 <- lme4breeding::lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_dam),
  relmat  = list(scanID = myKI),
  verbose = TRUE,
  data    = mydat_ordered
)
utils::summary(mix2)
vc <- lme4breeding::VarCorr(mix2)
print(vc, comp = c("Variance"))
ve <- attr(lme4breeding::VarCorr(mix2), "sc")^2
h2 <- vc$scanID / (vc$scanID + vc$F1_dam + ve)
base::as.numeric(h2)
vc$F1_dam / (vc$scanID + vc$F1_dam + ve)

mix3 <- lme4breeding::lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_sire) + (1|F1_dam),
  relmat  = list(scanID = myKI),
  verbose = TRUE,
  data    = mydat_ordered
)
utils::summary(mix3)
vc <- lme4breeding::VarCorr(mix3)
print(vc, comp = c("Variance"))
ve <- attr(lme4breeding::VarCorr(mix3), "sc")^2
h2 <- vc$scanID / (vc$scanID + vc$F1_sire + vc$F1_dam + ve)
base::as.numeric(h2)
vc$F1_dam / (vc$scanID + vc$F1_sire + vc$F1_dam + ve)
vc$F1_sire / (vc$scanID + vc$F1_sire + vc$F1_dam + ve)

# ---- dam 임의효과 행렬(K_dam) 구성 -------------------------------------
utils::str(mydat_ordered)
dam_ids <- base::factor(mydat_ordered$F1_dam)
Z <- stats::model.matrix(~ dam_ids - 1)
rownames(Z) <- mydat_ordered$F2_id

K_dam <- base::tcrossprod(Z)
K_dam <- K_dam + base::diag(1e-4, nrow(K_dam), ncol(K_dam))
base::dim(K_dam)

base::all.equal(rownames(K_dam), colnames(K_dam))
base::all.equal(colnames(K_dam), colnames(X_bw8_subset))
base::all.equal(mydat_ordered$F2_id, colnames(K_dam))
rownames(K_dam) <- mydat_ordered$scanID
colnames(K_dam) <- mydat_ordered$scanID

base::all.equal(rownames(K_dam), rownames(myKI))
base::all.equal(colnames(K_dam), colnames(myKI))
utils::str(bim_bw8_corrected)
base::all.equal(base::as.character(bim_bw8_corrected$index), rownames(X_bw8_subset))

# ---- Step 6: 널(null) 모형 적합 ----------------------------------------
columnIndex <- base::match("F2_BW8", colnames(scanAnno))
colnames(scanAnno)[columnIndex]

nullmod <- GENESIS::fitNullModel(                                         # [GENESIS::fitNullModel]
  scanAnno,
  outcome = colnames(scanAnno)[columnIndex],
  covars  = "F2_SEX",
  cov.mat = list(gen = myKI, dam = K_dam),
  family  = "gaussian"
)
nullmod$varComp
GENESIS::varCompCI(nullmod, prop = TRUE)                                  # [GENESIS::varCompCI]

nullmod2 <- GENESIS::fitNullModel(
  scanAnno,
  outcome = colnames(scanAnno)[columnIndex],
  covars  = c("F2_SEX","F1_dam"),
  cov.mat = myKI,
  family  = "gaussian"
)

# ---- Step 5): GWAS 연관성 검정 -----------------------------------------
genoIterator <- GWASTools::GenotypeBlockIterator(genoData, snpBlock = 10000)   # [GWASTools::GenotypeBlockIterator]

assoc_K_dam <- GENESIS::assocTestSingle(                                  # [GENESIS::assocTestSingle]
  genoIterator,
  null.model = nullmod,
  BPPARAM    = BiocParallel::SerialParam()                                # [BiocParallel::SerialParam]
)

assoc_fixed_dam <- GENESIS::assocTestSingle(
  genoIterator,
  null.model = nullmod2,
  BPPARAM    = BiocParallel::SerialParam()
)

base::save.image("ogye_bw8_gwas_updated.RData")

# ---- Step 8: P-값 보정 & 유의 SNP 선별 ---------------------------------
assoc_K_dam$p_adj.BH        <- stats::p.adjust(assoc_K_dam$Score.pval, method="BH")
assoc_K_dam$p_adj.bonferroni<- stats::p.adjust(assoc_K_dam$Score.pval, method="bonferroni")
assoc_sig_K_dam <- dplyr::filter(assoc_K_dam, p_adj.BH < 0.05)
assoc_sig_K_dam <- dplyr::left_join(assoc_sig_K_dam, bim_bw8_corrected, by = c("variant.id"="index"))
assoc_sig_K_dam_vep <- dplyr::left_join(assoc_sig_K_dam, bim_bw8_annotated, by = "id")

assoc_fixed_dam$p_adj.BH         <- stats::p.adjust(assoc_fixed_dam$Score.pval, method="BH")
assoc_fixed_dam$p_adj.bonferroni <- stats::p.adjust(assoc_fixed_dam$Score.pval, method="bonferroni")
assoc_sig_fixed_dam <- dplyr::filter(assoc_fixed_dam, p_adj.BH < 0.05)
assoc_sig_fixed_dam <- dplyr::left_join(assoc_sig_fixed_dam, bim_bw8_corrected, by = c("variant.id"="index"))

# ---- Step 9: 외부 GCTA 결과와 비교(선택) -------------------------------
gcta_mlma_sig_SNPs       <- final_gwas_results_annotated |> dplyr::filter(gwas_trait=="bw8") |> dplyr::distinct(id)
assoc_sig_K_dam_SNPs     <- assoc_sig_K_dam |> dplyr::distinct(id)
common_snps              <- base::intersect(gcta_mlma_sig_SNPs$id, assoc_sig_K_dam_SNPs$id)
assoc_sig_fixed_dam_SNPs <- assoc_sig_fixed_dam |> dplyr::distinct(id)
common_snps2             <- base::intersect(gcta_mlma_sig_SNPs$id, assoc_sig_fixed_dam_SNPs$id)

utils::write.csv(assoc_sig_K_dam, "bw8_gwas_results.csv", row.names = FALSE)

# ---- Step 10: QQ 플롯 & 원형 맨해튼 플롯 -------------------------------
library(CMplot)                                                            # [CMplot]

pvals  <- assoc_K_dam$Score.pval
chisq  <- stats::qchisq(1 - pvals, df = 1)
lambda <- stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, 1)
qqman::qq(pvals,                                                          # [qqman::qq]
          main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
          xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")

grDevices::png("bw8_gwas_qqplot.png", width=6, height=6, units="in", res=300)
qqman::qq(pvals, main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"))
grDevices::dev.off()

manplot <- assoc_K_dam |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)
base::names(manplot)

CMplot::CMplot(                                                           # [CMplot::CMplot]
  Pmap = manplot,
  plot.type = "q",
  conf.int = TRUE,
  box = TRUE,
  main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
  file.output = TRUE,
  file = "tiff"
)

manplot <- assoc_K_dam |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)

manplot <- dplyr::mutate(manplot, Chromosome = dplyr::case_when(Chromosome == "X" ~ "Z", TRUE ~ Chromosome))
dplyr::distinct(manplot, Chromosome)

raw_p <- manplot$P
adj_p <- stats::p.adjust(raw_p, method = "BH")
bh_threshold <- ifelse(any(adj_p <= 0.05), max(raw_p[adj_p <= 0.05], na.rm = TRUE), 0)
base::print(bh_threshold)

CMplot::CMplot(
  manplot,
  plot.type = "c",
  LOG10 = TRUE,
  threshold = c(bh_threshold, 0.05 / nrow(manplot)),
  threshold.col = c("orange","red"),
  chr.den.col = NULL,
  file = "jpg",
  dpi = 300, width = 9, height = 9
)

CMplot::CMplot(
  Pmap = manplot,
  plot.type = "d",
  bin.size = 1e6,
  file.output = TRUE,
  file = "jpg"
)

# --- (선택) 다형질 시각화 예시 -----------------------------------------
manplot$Chromosome <- base::as.character(manplot$Chromosome)
manplot$Chromosome[manplot$Chromosome == "M"] <- "34"
manplot$Chromosome[manplot$Chromosome == "X"] <- "36"
manplot$Chromosome <- base::as.numeric(manplot$Chromosome)
utils::str(manplot$Chromosome)

# 주의: multi_manhattan()는 본 스크립트 어디에서도 패키지를 로드하지 않았습니다.
# 보통 qqman/CMplot에는 이 함수가 없으니, 사용자 정의 함수이거나 별도 패키지일 수 있습니다.
# → 출처 확인 필요! 아래 라인은 출처 확인 후에만 사용하세요.
# multi_manhattan( ... )  # [출처 미상: 사용자 정의 함수 또는 별도 패키지]

# ---- 유의 SNP 유전자 주석 ------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"), ask = FALSE)
install.packages("dplyr")

library(rtracklayer); library(GenomicRanges); library(dplyr)

gtf_file <- "../Gallus_gallus_gca000002315v5.GRCg6a.114.gtf.gz"
gtf  <- rtracklayer::import(gtf_file)                                     # [rtracklayer::import]
genes<- gtf[gtf$type == "gene"]
utils::head(genes)

window_size <- 10000
assoc_sig_K_dam <- assoc_sig_K_dam %>%
  dplyr::mutate(chr = base::as.character(chr.x)) %>%
  dplyr::mutate(chr = base::ifelse(chr.x == "X","Z", chr.x))
assoc_sig_K_dam <- assoc_sig_K_dam %>%
  dplyr::mutate(chr = base::as.character(chr)) %>%
  dplyr::mutate(chr = base::ifelse(chr == "M","MT", chr))
dplyr::distinct(assoc_sig_K_dam, chr)

snp_ranges <- GenomicRanges::makeGRangesFromDataFrame(                    # [GenomicRanges::makeGRangesFromDataFrame]
  assoc_sig_K_dam %>%
    dplyr::mutate(start = pos.x - window_size,
                  end   = pos.x + window_size,
                  strand = "*"),
  keep.extra.columns = TRUE,
  na.rm = TRUE,
  seqnames.field = "chr"
)

overlaps <- GenomicRanges::findOverlaps(snp_ranges, genes, type = "any")  # [GenomicRanges::findOverlaps]
gene_summary <- data.frame(
  snp_index  = S4Vectors::queryHits(overlaps),                            # [S4Vectors::queryHits]
  gene_index = S4Vectors::subjectHits(overlaps)                           # [S4Vectors::subjectHits]
) %>%
  dplyr::mutate(
    snp_id      = assoc_sig_K_dam$id[snp_index],
    gene_id     = BiocGenerics::mcols(genes)$gene_id[gene_index],         # [BiocGenerics::mcols]
    gene_name   = BiocGenerics::mcols(genes)$gene_name[gene_index],
    gene_biotype= BiocGenerics::mcols(genes)$gene_biotype[gene_index],
    gene_start  = GenomicRanges::start(genes)[gene_index],                # [GenomicRanges::start]
    gene_end    = GenomicRanges::end(genes)[gene_index],                  # [GenomicRanges::end]
    gene_strand = base::as.character(GenomicRanges::strand(genes)[gene_index]) # [GenomicRanges::strand]
  ) %>%
  dplyr::group_by(snp_id) %>%
  dplyr::summarise(
    genes   = paste(unique(gene_name), collapse = "; "),
    details = paste(gene_id, " (", gene_start, "-", gene_end,
                    ", strand: ", gene_strand, ", biotype: ", gene_biotype, ")", collapse = "; ")
  )

annotated_snps <- assoc_sig_K_dam %>%
  dplyr::left_join(gene_summary, by = c("id" = "snp_id")) %>%
  dplyr::mutate(genes = base::ifelse(is.na(genes), "윈도우 내 유전자 없음", genes))
utils::head(annotated_snps)
utils::write.csv(annotated_snps, "annotated_significant_snps_bw8_window1Mb.csv", row.names = FALSE)

# ---- QC / rsID 주석 -----------------------------------------------------
library(dplyr); library(stringr)
raw_data     <- genio::read_plink("merge_total_chicken")                  # [genio::read_plink]
bim_raw_data <- raw_data$bim
dplyr::distinct(bim_raw_data, chr)
utils::head(bim_raw_data)

# chr=0에 있는 MT 변이 발견 → 필터 필요(--not-chr 0)
bim_with_rs <- bim_raw_data %>%
  dplyr::mutate(rs_id = stringr::str_extract(id, "rs\\d+"))               # [stringr::str_extract]
base::sum(!is.na(bim_with_rs$rs_id))

utils::write.csv(bim_with_rs, "bim_with_rs_bw8.csv", row.names = FALSE)

bim_with_rs %>%
  dplyr::filter(!is.na(rs_id)) %>%
  dplyr::select(rs_id) %>%
  utils::write.table("bim_bw8_rs_ids.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

vep_results_bim_with_rs_bw8 <- utils::read.delim("bim_with_rs_ids_vep.annotated.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(vep_results_bim_with_rs_bw8)[1] <- "rs_id"

bim_bw8_annotated <- dplyr::left_join(bim_with_rs, vep_results_bim_with_rs_bw8, by = "rs_id")
dplyr::distinct(bim_bw8_annotated, chr) %>% print(n=34)
bim_bw8_annotated |> dplyr::filter(chr == "LGE22C19W28_E50C23") |> dplyr::distinct(Location) %>% print(n=100)

# VEP 주석을 근거로 chr ID 보정
id1 <- bim_bw8_corrected |> dplyr::filter(chr == "LGE22C19W28_E50C23") |> dplyr::distinct(id)

bim_bw8_corrected <- bim_bw8_corrected %>%
  dplyr::mutate(chr = dplyr::case_when(
    chr == "Z" ~ "34",
    chr == "LGE22C19W28_E50C23" ~ "33",
    TRUE ~ chr
  ))

id2 <- bim_bw8_corrected |> dplyr::filter(chr == "33") |> dplyr::distinct(id)
base::all.equal(id1, id2)
dplyr::distinct(bim_bw8_corrected, chr)

utils::write.csv(assoc_sig_K_dam_vep, "bw8_gwas_results_updated.csv", row.names = FALSE)
