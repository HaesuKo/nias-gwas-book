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
# 4) GWASTools를 이용해 GWAS용 데이터 구조를 준비한다.
# 5) 개체 간 유전적 관계/집단 구조(GRM 기반 임의효과), 고정효과, 공변량을 반영한 귀무 모형(null model)을 적합한다.
# 6) 귀무 모형을 사용해 SNP 유전력(heritability)을 추정한다.
# 7) 유전자형과 표현형 간 연관성 검정(GWAS)을 수행한다.
# 8) p-값을 보정하고 유의한 SNP를 선별한다.

# ---- STEP 1. 라이브러리 로드 -------------------------------------------------
library(genio)         # PLINK I/O                [genio]
library(GENESIS)       # 귀무모형/연관성검정        [GENESIS]
library(GWASTools)     # 유전자형 데이터구조/iterator [GWASTools]
library(qqman)         # QQ/맨해튼 플롯           [qqman]
library(CMplot)         # QQ/맨해튼 플롯           [CMplot]
library(lme4breeding)  # LMM, 임의효과 표현형 분산비율 [lme4breeding]
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
utils::str(mydat)     
# [utils::str]
common_samples  <- base::intersect(colnames(X_bw8_corrected), mydat$F2_id)
X_bw8_subset    <- X_bw8_corrected[, common_samples]
mydat_ordered   <- mydat[base::match(common_samples, mydat$F2_id), ]
          
base::all.equal(mydat_ordered$F2_id, colnames(X_bw8_subset)) 
mydat_ordered$scanID <-  rownames(mydat_ordered)

# ---- Step 4: GRM 계산 --------------------------------------------------

geno        <- base::t(X_bw8_subset)
geno_coded  <- geno - 1
utils::str(geno_coded)

myKI <- sommer::A.mat(geno_coded)                                         # [sommer::A.mat]
myKI <- myKI + base::diag(1e-4, ncol(myKI), ncol(myKI))                   # [base::diag]

myKI[1:5, 1:5]
# Visualize GRM
colfunc <- grDevices::colorRampPalette(c("steelblue4","springgreen","yellow")) # [grDevices::colorRampPalette]
grDevices::png("bw8_Amatrix_updated.png", width=10, height=6, units="in", res=300) # [grDevices::png]
stats::heatmap(myKI, col = colfunc(100), Colv = "Rowv", symmetric = TRUE) # [stats::heatmap]
grDevices::dev.off()                                                      # [grDevices::dev.off]

# ---- Step 5: 데이터 간 변수 정렬 확인 & LMM 적합 ----------------------------------

colnames(myKI) = mydat_ordered$scanID
rownames(myKI) = mydat_ordered$scanID

# Check variable type
str(mydat_ordered)

# Dam as random effect
mix <- lme4breeding::lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_dam),
  relmat  = list(scanID = myKI),
  verbose = TRUE,
  data    = mydat_ordered
)

summary(mix)
vc <- VarCorr(mix)
print(vc, comp = c("Variance"))
ve <- attr(VarCorr(mix), "sc")^2
h2 <- vc$scanID / (vc$scanID + vc$F1_dam + ve) # SNP heritability 계산
base::as.numeric(h2)
dam_var = vc$F1_dam / (vc$scanID + vc$F1_dam + ve) #전체 표현형 분산 중 모계효과(Dam ID) 분산 비율 계산
base::as.numeric(dam_var)

# ---- dam 임의효과 행렬(K_dam) 구성 -------------------------------------
str(mydat_ordered)
dam_ids <- base::factor(mydat_ordered$F1_dam)
Z <- stats::model.matrix(~ dam_ids - 1)
rownames(Z) <- mydat_ordered$scanID

K_dam <- base::tcrossprod(Z)
K_dam <- K_dam + base::diag(1e-4, nrow(K_dam), ncol(K_dam))
base::dim(K_dam)


colnames(X_bw8_subset) = mydat_ordered$scanID

base::all.equal(colnames(K_dam), colnames(X_bw8_subset))
base::all.equal(rownames(K_dam), rownames(myKI))
base::all.equal(colnames(K_dam), colnames(myKI))


# ---- Step 6: GWASTools 데이터 구조 준비 --------------------------------
bim_bw8_corrected <- base::cbind(index = 1:nrow(bim_bw8_corrected), bim_bw8_corrected)
str(bim_bw8_corrected)

base::all.equal(rownames(X_bw8_subset), bim_bw8_corrected$id)

rownames(X_bw8_subset) <- bim_bw8_corrected$index

all.equal(as.character(bim_bw8_corrected$index), rownames(X_bw8_subset))

dplyr::distinct(bim_bw8_corrected, chr) %>% print()

rownames(X_bw8_subset)
colnames(X_bw8_subset)
all.equal(colnames(X_bw8_subset), mydat_ordered$scanID)

# 종에 따라 염색체 코딩을 다르게 해야한다:autosomeCode, XchromCode, etc.
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

# ---- Step 7: 귀무 모형 적합 및 유전력 계산 ---------------------------
columnIndex <- base::match("F2_BW8", colnames(scanAnno))
colnames(scanAnno)[columnIndex]

nullmod <- GENESIS::fitNullModel(                                         # [GENESIS::fitNullModel]
  scanAnno,
  outcome = colnames(scanAnno)[columnIndex],
  covars  = "F2_SEX",
  cov.mat = list(gen = myKI, dam = K_dam),
  family  = "gaussian"
)

# SNP heritability 계산 (h2, CI value)
nullmod$varComp
GENESIS::varCompCI(nullmod, prop = TRUE)                                  # [GENESIS::varCompCI]

# ---- Step 8: GWAS 연관성 검정 -----------------------------------------
genoIterator <- GWASTools::GenotypeBlockIterator(genoData, snpBlock = 10000)   # [GWASTools::GenotypeBlockIterator]

assoc_K_dam <- GENESIS::assocTestSingle(                                  # [GENESIS::assocTestSingle]
  genoIterator,
  null.model = nullmod,
  BPPARAM    = BiocParallel::SerialParam()                                # [BiocParallel::SerialParam]
)

base::save.image("ogye_bw8_gwas_updated.RData")

# ---- Step 9: P-값 보정 & 유의 SNP 선별 ---------------------------------
assoc_K_dam$p_adj.BH        <- stats::p.adjust(assoc_K_dam$Score.pval, method="BH")
assoc_K_dam$p_adj.bonferroni<- stats::p.adjust(assoc_K_dam$Score.pval, method="bonferroni")
assoc_sig_K_dam <- dplyr::filter(assoc_K_dam, p_adj.BH < 0.05)
assoc_sig_K_dam <- dplyr::left_join(assoc_sig_K_dam, bim_bw8_corrected, by = c("variant.id"="index"))
assoc_sig_K_dam_vep <- dplyr::left_join(assoc_sig_K_dam, bim_bw8_annotated, by = "id")

# ---- Step 10: QQ 플롯 & 원형 맨해튼 플롯 -------------------------------
library(CMplot)                                                            # [CMplot]

pvals  <- assoc_K_dam$Score.pval
chisq  <- stats::qchisq(1 - pvals, df = 1)
lambda <- stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, 1)

manplot <- assoc_K_dam |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)
base::names(manplot)

manplot <- dplyr::mutate(manplot, Chromosome = dplyr::case_when(Chromosome == "X" ~ "Z", TRUE ~ Chromosome))
dplyr::distinct(manplot, Chromosome)

## QQ plot
CMplot::CMplot(                                                           # [CMplot::CMplot]
  Pmap = manplot,
  plot.type = "q",
  conf.int = TRUE,
  box = TRUE,
  main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
  file.output = TRUE,
  file = "tiff"
)

## Circular Manhattan plot
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

## SNP density plot across chromosomes
CMplot::CMplot(
  Pmap = manplot,
  plot.type = "d",
  bin.size = 1e6,
  file.output = TRUE,
  file = "jpg"
)


###

rmarkdown::render("gwas_workflow_vF.Rmd")
browseURL("gwas_workflow_vF.html")

rmarkdown::render("2.gwas_workflow.Rmd")
browseURL("2.gwas_workflow.html")

### 
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata#[,1:200]
MP <- MP_cpdata
M<- GT
n <- nrow(DT) # to be used for degrees of freedom
k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)


```{=html}
<style>
  /* ===== GWAS readable UI (conflict-safe) ===== */
  
  /* 카드/컬러 — 변수 없이 직접 지정(충돌 회피) */
  .info-box{
    padding: 1em 1.2em;
    margin: 1.1em 0;
    background-color: #e7f3ff !important;
      border-left: 5px solid #4a90e2 !important;
    border-radius: 8px;
    color: #15457a !important;
  }
.info-box h4, .info-box h5{
  margin: 0 0 .4rem;
  color: #4a90e2 !important;
}
.info-box.compact{ font-size: .95rem; line-height: 1.5; }

.info-box.success{
  background-color: #e6f4ea !important;
    border-left-color: #34a853 !important;
    color: #0b3d02 !important;
}
.info-box.success h4, .info-box.success h5{ color: #34a853 !important; }
    
    .info-box.warn{
      background-color: #fff7e6 !important;
        border-left-color: #f59e0b !important;
        color: #7c3e00 !important;
    }
  .info-box.warn h4, .info-box.warn h5{ color: #f59e0b !important; }
      
      /* 그리드 카드 */
      .gw-grid{
        display: grid !important;
        gap: 12px !important;
        margin: 12px 0 !important;
        grid-template-columns: 1fr !important;
      }
    @media(min-width:900px){
      .gw-grid{ grid-template-columns: repeat(2, 1fr) !important; }
    }
    .gw-card{
      background: #fff !important;
        border: 1px solid #e5e7eb !important;
      border-radius: 10px !important;
      padding: .9rem 1rem !important;
      color: #374151 !important;
        box-shadow: 0 1px 2px rgba(0,0,0,.04) !important;
    }
    .gw-card h4, .gw-card h5{ margin:.1rem 0 .5rem !important; }
    
    /* 접기 패널(details) */
      .gw-toggle{
        border:1px solid #d1d5db !important;
        border-radius:10px !important;
        background:#fff !important;
          margin:.7rem 0 !important;
        overflow:hidden !important;
      }
    .gw-toggle > summary{
      cursor:pointer !important;
      list-style:none !important;
      padding:.8rem 1rem !important;
      font-weight:700 !important;
      background:#f9fafb !important;
        color:#111827 !important;
        display:flex !important;
      align-items:center !important;
      justify-content:space-between !important;
    }
    .gw-toggle > summary::marker,
    .gw-toggle > summary::-webkit-details-marker{ display:none !important; }
    .gw-toggle[open] > summary{ border-bottom:1px solid #e5e7eb !important; }
      .gw-toggle .gw-toggle-body{ padding:.85rem 1rem !important; color:#374151 !important; }
          
          /* 컴팩트 리스트/주석 */
          .compact-list{ margin:.2rem 0 .1rem !important; padding-left:1.1rem !important; }
        .compact-list li{ margin:.15rem 0 !important; }
        .small-note{ font-size:.9rem !important; color:#4b5563 !important; }
            
            /* 코드칩 톤 맞춤 */
            .info-box code, .gw-card code, .gw-toggle-body pre{
              background:#f8f9fa !important; padding:.08rem .28rem !important; border-radius:4px !important;
            }
          .gw-toggle-body pre { padding: .5rem .7rem !important; }
          </style>
            
            <h3>lme4 vs. lme4breeding: 핵심 차이점 비교 🧬</h3>
            
            <div class="gw-grid">
            <div class="gw-card">
            <h4>1) lme4::lmer (범용 혼합 모델)</h4>
  <p class="small-note">일반적인 다수준/반복측정 데이터 분석용</p>
    <div class="info-box compact">
      <h5>핵심 가정</h5>
      <p>임의 효과로 지정된 개체(또는 그룹)들은 기본적으로 <strong>서로 독립적(independent)</strong>이라고 가정</p>
      </div>
      <ul class="compact-list">
        <li><strong>주요 목적</strong>: 그룹 구조에서 오는 분산을 통제</li>
        <li><strong>사용 분야</strong>: 사회과학, 의학, 생태학 등 광범위</li>
        <li><strong>한계</strong>: 개체 간 상관관계를 직접 지정하기 어려움</li>
        </ul>
        <details class="gw-toggle">
          <summary>기본 사용 예시</summary>
          <div class="gw-toggle-body">
            <p>학생(<code>Student</code>)과 학교(<code>School</code>) 효과를 고려하는 모델. 각 학생은 다른 학생과 독립적이라고 가정.</p>
            <pre><code># lme4 패키지 사용
            lmer(score ~ (1|Student) + (1|School), data = df)</code></pre>
            </div>
            </details>
            </div>
            
            <div class="gw-card">
              <h4>2) lme4breeding::lmebreed (육종 특화)</h4>
  <p class="small-note">유전적 관계를 고려한 분석용</p>
    <div class="info-box success compact">
      <h5>핵심 기능</h5>
      <p><code>relmat</code> 인자로 <strong>유전 관계 행렬(Kinship Matrix)</strong>을 모델에 직접 반영하여, 혈연관계를 통계적으로 고려</p>
      </div>
      <ul class="compact-list">
        <li><strong>주요 목적</strong>: 유전력 및 육종가(breeding value) 추정</li>
        <li><strong>사용 분야</strong>: 동물/식물 육종, 유전체학</li>
        <li><strong>강점</strong>: 집단 구조 효과를 정교하게 보정 가능</li>
        </ul>
        <details class="gw-toggle">
          <summary>혈연관계 적용 예시</summary>
          <div class="gw-toggle-body">
            <p>개체(<code>scanID</code>)의 유전적 관계(<code>K_scan</code>)를 반영하여 유전 효과를 추정. 친척은 비슷한 효과를 공유한다고 가정.</p>
            <pre><code># lme4breeding 패키지 사용
            lmebreed(phenotype ~ (1|scanID), 
                     relmat = list(scanID = K_scan), 
                     data = PH)</code></pre>
            </div>
            </details>
            </div>
            </div>
            ```
          