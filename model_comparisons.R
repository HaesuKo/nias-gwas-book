## -------------------------
## Model comparison template
## -------------------------
suppressPackageStartupMessages({
  library(lme4breeding)
  library(Matrix)
})

## ---- 0) Prep helpers ----------------------------------------------------
drop_na_outcome <- function(df, y) {
  df <- df[!is.na(df[[y]]), , drop = FALSE]
  droplevels(df)
}

align_kernel_levels <- function(K, lev) {
  stopifnot(setequal(rownames(K), colnames(K)))
  K[lev, lev, drop = FALSE]
}

ridge <- function(K, eps = 1e-4) {
  diag(K) <- diag(K) + eps
  K
}

summ_vc <- function(mod) {
  vc <- VarCorr(mod)
  ve <- attr(vc, "sc")^2
  list(vc = vc, ve = ve)
}

safe_aicbic <- function(mod) {
  c(AIC = AIC(mod), BIC = BIC(mod), logLik = as.numeric(logLik(mod)))
}

## ---- 1) Data & kernels ---------------------------------------------------
## REQUIRED objects in your workspace:
##   pheno_holstein, myKI
## pheno_holstein has columns days_in_milk, scanID, farm_id.
## myKI is a kinship (scanID×scanID) with row/colnames covering the used scanIDs.
## OPTIONAL:
##   K_farmL  (farm-by-farm covariance)

YVAR <- "days_in_milk"
PH <- droplevels(subset(pheno_holstein, !is.na(days_in_milk)))
PH$scanID  <- factor(PH$scanID)
PH$farm_id <- factor(PH$farm_id)

K_scan <- myKI[levels(PH$scanID), levels(PH$scanID)]


## ---- 2) Fit models -------------------------------------------------------
## Pattern A: random farm (no farm kernel), genetic kernel for scanID

mix_A <- lme4breeding::lmebreed(
  days_in_milk ~ (1|scanID) + (1|farm_id),
  relmat  = list(scanID = K_scan),   # only for the genetic id effect
  data    = PH,
  verbose = TRUE
)


# farm kernel must be farm-by-farm
Zf <- model.matrix(~ 0 + PH$farm_id)           # n × L
lev_farm <- colnames(Zf)                        # farm levels
# EXAMPLE kernel: identity (change to what you need)
K_farmL <- diag(length(lev_farm))
dimnames(K_farmL) <- list(lev_farm, lev_farm)
# strip "PH$farm_id" from row/col names
nm <- rownames(K_farmL)
nm2 <- sub("^PH\\$farm_id", "", nm)
dimnames(K_farmL) <- list(nm2, nm2)

# align to the actual farm levels used in the data (and check)
lev <- levels(PH$farm_id)
stopifnot(identical(sort(nm2), sort(lev)))

# reorder to match the model’s factor level order (optional but nice)
K_farmL <- K_farmL[lev, lev, drop = FALSE]

mix_B <- lme4breeding::lmebreed(
  days_in_milk ~ (1|scanID) + (1|farm_id),
  relmat  = list(scanID = K_scan, farm_id = K_farmL),
  data    = PH,
  verbose = TRUE
)


## Observation-by-farm incidence (Pattern C)
Z_sameFarm <- model.matrix(~ 0 + farm_id, data = PH)  # n × L
colnames(Z_sameFarm) <- levels(PH$farm_id)            # clean
rownames(Z) <- as.character(scanID_int)
K_farm <- tcrossprod(Z)
K_farm <- K_farm + diag(1e-4, nrow(K_farm))
dimnames(K_farm) <- list(as.character(scanID_int), as.character(scanID_int))

## ---- 2) Fit models -------------------------------------------------------
## Pattern A: random farm (no farm kernel), genetic kernel for scanID
mix_A <- lmebreed(
  formula = as.formula(sprintf("%s ~ 1 + (1|scanID) + (1|farm_id)", YVAR)),
  relmat  = list(scanID = K_scan),
  data    = PH,
  verbose = FALSE
)

## Pattern B: random farm WITH farm-level covariance K_farmL
mix_B <- lmebreed(
  formula = as.formula(sprintf("%s ~ 1 + (1|scanID) + (1|farm_id)", YVAR)),
  relmat  = list(scanID = K_scan, farm_id = K_farmL),
  data    = PH,
  verbose = FALSE
)

## Pattern C: custom incidence via addmat (farm_obs), genetic kernel for scanID
mix_C <- lmebreed(
  formula = as.formula(sprintf("%s ~ 1 + (1|scanID) + (1|farm_id)", YVAR)),
  relmat  = list(scanID = K_scan),
  addmat  = list(farm_obs = Z_sameFarm),
  data    = PH,
  verbose = FALSE
)

## Pattern : farm effect 200 X 200 
mix_D <- lmebreed(
  formula = as.formula(sprintf("%s ~ 1 + (1|scanID) + (1|farm_id)", YVAR)),
  relmat  = list(scanID = K_scan),
  addmat  = list(farm_id = Z_sameFarm),
  data    = PH,
  verbose = FALSE
)

## Farm as FIXED effect (+ genetic random intercept)
mix_fix <- lmebreed(
  formula = as.formula(sprintf("%s ~ farm_id + (1|scanID)", YVAR)),
  relmat  = list(scanID = K_scan),
  data    = PH,
  verbose = FALSE
)

## ---- 3) Collect fit stats ------------------------------------------------
fits <- list(
  Pattern_A = mix_A,
  Pattern_B = mix_B,
  Pattern_C = mix_C,
  pattern_D = mix_D,
  FarmFixed = mix_fix
)

fit_tab <- do.call(rbind, lapply(fits, safe_aicbic))
fit_tab <- as.data.frame(fit_tab)
fit_tab$k  <- sapply(fits, function(m) attr(logLik(m), "df"))
fit_tab$n  <- nrow(PH)
fit_tab$REML <- TRUE  # lmebreed uses REML by default

## ---- 4) Variance components snapshot ------------------------------------
vc_summ <- lapply(fits, function(m) {
  s <- summ_vc(m)
  list(
    Ve = unname(s$ve),
    Rand = s$vc
  )
})

## ---- 5) Print results ----------------------------------------------------
cat("\n=== Model Fit Comparison (REML) ===\n")
print(round(fit_tab, 3))

cat("\nNote:\n- AIC/BIC are comparable across A/B/C (same fixed effects, REML).\n",
    "- Comparing FarmFixed versus others is safer with ML; if desired, refit with REML=FALSE.\n")

cat("\n=== Residual Variance (Ve) by model ===\n")
print(sapply(vc_summ, function(x) x$Ve))

cat("\n=== Random-effect variance components (glance) ===\n")
for (nm in names(vc_summ)) {
  cat(sprintf("\n-- %s --\n", nm))
  print(vc_summ[[nm]]$Rand, comp = "Variance")
}

## ---- 6) (Optional) ML refits for fixed-effect comparison -----------------
## If you want to compare FarmFixed vs Pattern_A via likelihood ratio:
## mix_A_ML   <- update(mix_A, REML = FALSE)
## mix_fix_ML <- update(mix_fix, REML = FALSE)
## anova(mix_A_ML, mix_fix_ML)  # LRT (different fixed effects -> use ML!)
