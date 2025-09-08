## 1. 패키지 설치 & 로딩
# ---- 설치어시스턴트 -------------------------------------------------------
cran_pkgs <- c("genio", "qqman", "sommer", "lme4breeding", "dplyr", "CMplot", "stringr", "ggplot2", "tibble", "readr")
bioc_pkgs <- c("GENESIS","GWASTools","SNPRelate","BiocParallel","rtracklayer","GenomicRanges")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install_if_missing <- function(pkgs, installer) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) installer(need, ask = FALSE)
}

install_if_missing(cran_pkgs, install.packages)
install_if_missing(bioc_pkgs, BiocManager::install)

# ---- 라이브러리 로드 ------------------------------------------------------
library(genio)
library(GENESIS)
library(GWASTools)
library(qqman)
library(CMplot)
library(lme4breeding)
library(sommer)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)

## 2. 경로/입력 정의
breed        <- "holstein"
plink_prefix <- "03_qc/NIAS_ibv3_296ea.holstein"
pheno_file   <- "NIAS_ibv3_296ea_pheno.csv"

report_dir <- file.path("05_gwas", breed)
figs_dir   <- file.path(report_dir, "figs")
gwas_dir   <- file.path(report_dir, "gwas")
tables_dir <- file.path(report_dir, "tables")
for (d in c(report_dir, figs_dir, gwas_dir, tables_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## 3. PLINK binary 데이터 불러오기
# read_plink()는 SNP x Sample 행렬(X), bim(마커정보), fam(개체정보)을 불러온다.
plink <- read_plink(plink_prefix)

X <- plink$X             # 행: SNP, 열: 샘플, 값: 0/1/2
bim <- plink$bim         # chr, id, pos, alt, ref 등
fam <- plink$fam         # family/sample id 등

# 미리보기
X[1:5, 1:5]
head(bim)
head(fam)

## 4. 표현형 데이터 로딩 & 정렬
pheno <- readr::read_csv(pheno_file, show_col_types = FALSE)

sample_id_chr <- colnames(X)  # PLINK 열이 샘플 ID
scanID <- seq_along(sample_id_chr)

pheno_ord <- pheno[match(sample_id_chr, pheno$individual_id), , drop = FALSE]
stopifnot(!any(is.na(pheno_ord$individual_id)))
pheno_ord$sample_id <- sample_id_chr
pheno_ord$scanID    <- scanID
stopifnot(identical(as.character(pheno_ord$individual_id), pheno_ord$sample_id))
# 미리보기
str(pheno_ord)

# 표현형 분포 확인
# --- 분석할 형질을 이 변수에 지정 ---
single_trait <- "days_in_milk"

p=ggplot(pheno_ord, aes(x = .data[[single_trait]])) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey70", color = "white") +
  geom_density(linewidth = 1) +
  labs(
    title = paste("Distribution of", single_trait),
    x = single_trait, 
    y = "Density"
  ) +
  theme_minimal()

# 변수에 저장된 플롯 확인
print(p)
# ggsave() 함수를 사용하여 플롯을 png 파일로 저장
ggsave(
  filename = file.path(figs_dir, paste0("pheno_distribution_", single_trait, ".png")),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)
## 5. GRM(유전체 관계 행렬) 계산 & 간단 시각화
# A.mat은 -1/0/1 코딩된 genotype matrix 이용
G_mat <- t(X)              # 샘플 x SNP
colnames(G_mat) <- rownames(X)
rownames(G_mat) <- sample_id_chr

G_code <- G_mat - 1        # 0/1/2 -> -1/0/1
KI     <- sommer::A.mat(G_code)
# 수치안정성(아주 작은 ridge)
diag(KI) <- diag(KI) + 1e-4

# scanID 기준의 dimnames 지정
rownames(KI) <- as.character(scanID)
colnames(KI) <- as.character(scanID)

stopifnot(
  identical(rownames(KI), as.character(scanID)),
  identical(colnames(KI), as.character(scanID))
)

# 상위 5x5 확인
KI[1:5, 1:5]

# 히트맵
colfunc <- grDevices::colorRampPalette(c("steelblue4","springgreen","yellow"))

png(file.path(report_dir, "A_matrix_heatmap.png"), width=2000, height=1400, res=200)
stats::heatmap(KI, col = colfunc(100), Colv = "Rowv", symmetric = TRUE)
dev.off()

## 6. GWAS 분석용 데이터 객체 만들기
bim <- cbind(index = seq_len(nrow(bim)), bim)
stopifnot(all.equal(rownames(X), bim$id))
stopifnot(all.equal(colnames(X), pheno_ord$sample_id))

rownames(X) <- bim$index
colnames(X) <- scanID

# GENESIS 패키지에서 요구하는 데이터 객체 만들기
mg <- GWASTools::MatrixGenotypeReader(
  genotype   = X,
  snpID      = bim$index,
  chromosome = as.integer(bim$chr),
  position   = as.integer(bim$pos),
  scanID     = scanID,
  autosomeCode = 1L:29L,
  XchromCode   = 30L, YchromCode = 31L, XYchromCode = 32L, MchromCode = 33L
)
# 유전자형 데이터 객체 생성
genoData <- GWASTools::GenotypeData(mg)
# 표현형 데이터 객체 생성
scanAnno <- GWASTools::ScanAnnotationDataFrame(pheno_ord)
# 객체 구조 확인
str(genoData)
str(scanAnno@data)   

## 7. GWAS 귀무모형 적합 및 유전력(SNP heritability) 계산
# --- 분석할 형질을 이 변수에 지정 ---
single_trait <- "days_in_milk"

# 분석할 형질이 데이터에 존재하는지 확인. 없으면 에러를 발생시키고 중단
stopifnot(single_trait %in% colnames(scanAnno@data))

# --- 데이터 전처리 ---

# 1. 모델에서 효과를 분석할 변수들을 factor(범주형 변수)로 변환
pheno_ord$farm_id <- as.factor(pheno_ord$farm_id)
pheno_ord$scanID <- as.factor(pheno_ord$scanID)

# 2. 분석할 형질(days_in_milk)에 결측치(NA)가 있는 샘플을 제거하여
#    최종 분석 데이터셋 'PH'를 생성
#    droplevels() 함수는 PH 데이터셋에 더 이상 존재하지 않는 샘플이나 농장의 레벨 정보를 제거
PH <- droplevels(subset(pheno_ord, !is.na(days_in_milk)))

# 3. Kinship 행렬(KI)을 최종 분석 데이터셋 'PH'에 맞춰 자른다.
#    반드시 'PH'의 샘플 목록과 순서를 기준으로 잘라야,
#    표현형 데이터와 Kinship 데이터의 샘플이 1:1로 정확하게 일치한다.
K_scan <- KI[levels(PH$scanID), levels(PH$scanID)]

# 농장 효과와 개체별 유전 효과를 임의 효과(random effect)로 설정하여 모델을 만든다.
# lmebreed 함수는 혈연관계(kinship)를 고려한 혼합 모델을 적합해준다.
mix <- lme4breeding::lmebreed(
  # 모델 공식: days_in_milk를 개체 효과(scanID)와 농장 효과(farm_id)로 설명
  days_in_milk ~ (1|scanID) + (1|farm_id),
  
  # relmat 인자: scanID의 임의 효과는 K_scan 행렬에 정의된 유전적 관계를 따르도록 지정한다.
  relmat  = list(scanID = K_scan),
  
  # verbose=TRUE: 모델 계산 과정을 화면에 출력
  verbose = TRUE,
  
  # data 인자: 모델에 사용할 데이터셋으로 PH를 지정
  data    = PH
)

# 모델 분석 결과 요약본 출력
summary(mix)

### 함수를 이용한 유전력 계산
# VarCorr() 함수로 모델('mix')에서 분산 성분(Variance Components)을 추출
vc <- VarCorr(mix)
print(vc, comp = c("Variance"))

# 잔차 분산(Residual variance)을 추출
# VarCorr 결과의 속성(attribute) "sc"는 잔차의 표준편차이므로 제곱하여 분산을 구한다. 
ve <- attr(VarCorr(mix), "sc")^2

# SNP 유전력(h2) 계산
# h2 = 유전분산 / (유전분산 + 농장분산 + 잔차분산)
h2 <- vc$scanID / (vc$scanID +  vc$farm_id + ve)
print(paste("SNP Heritability (h2):", base::as.numeric(h2)))

# 전체 표현형 분산 중 농장 효과(farm ID)가 차지하는 분산 비율 계산
farm_var <- vc$farm_id / (vc$scanID + vc$farm_id + ve)
print(paste("Farm Variance ratio:", base::as.numeric(farm_var)))

## 8. GWAS 최종 null 모형 적합 및 GWAS
# --- 분석할 형질을 이 변수에 지정 ---
single_trait <- "days_in_milk"

# 분석할 형질이 데이터에 존재하는지 확인. 없으면 에러를 발생시키고 중단
stopifnot(single_trait %in% colnames(scanAnno@data))

# 귀무모형
null_single <- GENESIS::fitNullModel(
  x        = scanAnno,
  outcome  = single_trait,
  cov.mat  = list(gen = KI),
  family   = "gaussian"
)

# 모형에서 추정된 분산 성분(유전분산, 잔차분산)을 확인
null_single$varComp 
# SNP heritability 계산 (h2, CI value)
GENESIS::varCompCI(null_single, prop = TRUE)       

# GWAS 실행
# 유전체 데이터를 **블록 단위(여기서는 10,000 SNP씩)**로 순차 처리하기 위한 이터레이터 생성: 메모리 사용을 줄이기 위함
genoIt <- GWASTools::GenotypeBlockIterator(genoData, snpBlock = 10000)
# 위에서 적합한 귀무모형을 고정한 상태에서, 각 SNP에 대해 단일변량 Score test로 연관성을 계산
assoc_single <- GENESIS::assocTestSingle(
  genoIt, null.model = null_single, BPPARAM = BiocParallel::SerialParam()
)

# ---- P-값 보정 & 유의 SNP 선별 ---------------------------------
assoc_single$p_adj_BH=p.adjust(assoc_single$Score.pval, method="BH")
assoc_single$p_adj_bonferroni=p.adjust(assoc_single$Score.pval, method="bonferroni")

# 저장
saveRDS(assoc_single, file = file.path(gwas_dir, paste0("assoc_", single_trait, ".rds")))
readr::write_csv(assoc_single, file.path(tables_dir, paste0("assoc_", single_trait, ".csv")))

head(assoc_single, 5)

# 결과 요약
# genomic inflation 구하기
chisq  <- stats::qchisq(1 - assoc_single$Score.pval, df = 1)
lambda <- median(chisq, na.rm = TRUE)/stats::qchisq(0.5, 1)
# 요약 표로 정리하기
summary_single <- tibble::tibble(
  trait = single_trait,
  n_snps = sum(!is.na(assoc_single$Score.pval)),
  lambda = lambda,
  n_sig_BH = sum(assoc_single$p_adj_BH < 0.05, na.rm = TRUE),
  n_sig_bonf = sum(assoc_single$p_adj_bonferroni < 0.05, na.rm = TRUE),
  min_p = min(assoc_single$Score.pval, na.rm = TRUE),
  min_p_SNP = assoc_single$variant.id[which.min(assoc_single$Score.pval)],
  min_BH = min(assoc_single$p_adj_BH, na.rm = TRUE)
)
readr::write_csv(summary_single, file.path(tables_dir, paste0("gwas_summary_", single_trait, ".csv")))

## 9. GWAS 결과 시각화
# PVE 분포 확인
ggplot(assoc_single, aes(x = PVE*100)) + # PVE를 퍼센트로 변환
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey70", color = "white") +
  geom_density(linewidth = 1) +
  labs(title = "Distribution of Variance Explained (PVE) by Single Variants", x = "PVE (%)", y = "Density") +
  theme_minimal()

# 플롯 저장
ggsave(file.path(figs_dir, "pve_hist.png"), dpi=300, width=6, height=4)

# ---- QQ 플롯 & 맨해튼 플롯 -------------------------------
manplot <- assoc_single |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)

# 염색체 분류 확인
dplyr::distinct(manplot, Chromosome)

# U 염색체 코딩 변경
manplot <- dplyr::mutate(manplot, Chromosome = dplyr::case_when(Chromosome == "U" ~ "X", TRUE ~ Chromosome))

## QQ plot
CMplot::CMplot(                                                           #
  Pmap = manplot,
  plot.type = "q",
  conf.int = TRUE,
  box = TRUE,
  main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
  file = "png",
  dpi = 300, width = 4, height = 4,
  # ↓↓↓ 크기 조절
  axis.cex = 0.7,   # 축 눈금 글자 크기(기본 1 → 더 작게)
  lab.cex  = 0.9,   # 축 라벨(제목) 크기(기본 1.5 → 더 작게)
  main.cex = 1.0,   # 타이틀 크기(기본 1.5 → 더 작게)
  # 라벨 위치 미세조정
  xticks.pos = 0.5,    # x축 눈금과 축 간격
  ylab.pos   = 1.5,    # y축 라벨과 축 간격
  mar = c(3, 4, 2.5, 2)  # 아래-왼-위-오른 여백(클립 방지/밀착 조절)
)

## Manhattan plot
CMplot::CMplot(
  manplot,
  plot.type = "m",
  LOG10 = TRUE,
  chr.den.col = NULL,
  file = "png",
  dpi = 300, width = 12, height = 4
)

## Circular Manhattan plot
CMplot::CMplot(
  manplot,
  plot.type = "c",
  LOG10 = TRUE,
  chr.den.col = NULL,
  file = "png",
  dpi = 300, width = 4, height = 4,
  axis.cex = 0.7,   # 축 눈금 글자 크기(기본 1 → 더 작게)
  lab.cex  = 0.9,   # 축 라벨(제목) 크기(기본 1.5 → 더 작게)
  mar = c(1, 1, 1, 1),   # 기본(3,6,3,3)보다 훨씬 타이트
  cir.chr.h = 0.6              # 바깥 크로모솜 테두리 두께 축소
)
