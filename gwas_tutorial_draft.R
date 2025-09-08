# Step-by-Step Tutorial: GWAS Analysis in R for Body Weight Trait (BW8) Using Chicken Genotype Data

# This tutorial guides beginners through performing a Genome-Wide Association Study (GWAS) 
# in R using genotype data from PLINK files and phenotype data. We will:
# 1. Load necessary libraries and data.
# 2. Prepare genotype and phenotype data.
# 3. Compute the Genomic Relationship Matrix (GRM).
# 4. Prepare data structures for GWAS using GWASTools.
# 5. Fit a null model accounting for relatedness and covariates.
# 6. Run the association test.
# 7. Adjust p-values and identify significant SNPs.
# 8. Compare results with external (e.g., GCTA) results.
# 9. Estimate heritability using the GRM.

# Prerequisites:
# - Install required packages if not already done: install.packages(c("genio", "GENESIS", "SNPRelate", "GWASTools", "qqman", "sommer", "lme4breeding", "dplyr"))
# - Have PLINK files ('bw8_corrected.bed', 'bw8_corrected.bim', 'bw8_corrected.fam') ready.
# - Have phenotype data in a data frame called 'merged5' with columns like 'F2_id', 'F2_BW8', 'F2_SEX', 'F1_dam'.
# - Have external GWAS results in 'final_gwas_results_annotated' for comparison (optional).
# - Note: This assumes biallelic SNPs coded as 0/1/2 (allele counts) in the genotype matrix.
# - For chicken genome: Autosomes are chromosomes 1-33, with custom codes for sex chromosomes (Z=36, W=35, etc.).

# Step 1: Load Required Libraries
# These packages handle PLINK data (genio), GWAS tools (GENESIS, GWASTools), relationship matrices (sommer), 
# visualization (qqman), mixed models (lme4breeding), and data manipulation (dplyr).
library(genio)       # For reading PLINK files.
library(GENESIS)     # For fitting null models and association tests in GWAS.
library(GWASTools)   # For genotype data handling and iterators in GWAS.
library(qqman)       # For creating QQ and Manhattan plots (used implicitly later if needed).
library(lme4breeding) # For fitting linear mixed models for gwas (find covariates, random effects), and estimating heritability.
library(dplyr)       # For data manipulation like filtering and joining.

# Set working directory and check if files are present.
setwd("~/ogye_leghorn_gwas/example/")
list.files()

# Step 2: Load PLINK Files for Genotype Data
# Read the PLINK files into R objects: X (genotype matrix), bim (SNP info), fam (sample info).

bw8_input_corrected <- read_plink('bw8_corrected')  # Load PLINK data; assumes files are in working directory.
X_bw8_corrected <- bw8_input_corrected$X            # Extract genotype matrix (rows: SNPs, columns: samples, values: 0/1/2).
bim_bw8_corrected <- bw8_input_corrected$bim        # Extract BIM data frame (SNP details: chr, id, pos, alt, ref).
fam_bw8_corrected <- bw8_input_corrected$fam        # Extract FAM data frame (sample details).

# Step 3: Prepare Phenotype Data
# Assume 'merged5' is your full phenotype dataset with sample IDs in 'F2_id' and phenotype in 'F2_BW8'.
sample_ids <- colnames(X_bw8_corrected)  # Get sample IDs from genotype matrix columns.
mydat <- merged5[merged5$F2_id %in% sample_ids, ]  # Subset phenotypes to match genotyped samples.
head(mydat, 6)  # Preview the first few rows of phenotype data.
str(mydat)  # Check the structure of the phenotype data frame.

# Step 4: Compute the Genomic Relationship Matrix (GRM)
# Subset genotype matrix to match phenotype samples (resolve dimension mismatches).
common_samples <- intersect(colnames(X_bw8_corrected), mydat$F2_id)  # Find overlapping sample IDs.
X_bw8_subset <- X_bw8_corrected[, common_samples]  # Subset genotype matrix to common samples.
mydat_ordered <- mydat[match(common_samples, mydat$F2_id), ]  # Reorder phenotypes to match genotype order.
all.equal(mydat_ordered$F2_id, colnames(X_bw8_subset))
rownames(X_bw8_subset)

# Transpose and recode genotypes for sommer or lme4breeding's A.mat (expects rows: samples, columns: SNPs, coded -1/0/1).
geno <- t(X_bw8_subset)  # Transpose: now rows are samples, columns are SNPs.
head(geno)
geno_coded <- geno - 1  # Recode: 0 -> -1 (homozygous minor), 1 -> 0 (heterozygous), 2 -> 1 (homozygous major).
rownames(geno_coded)  # Extract sample IDs from rows.
# To make the structure of geno_coded identical to a plain numeric matrix without dimension names,
# remove the dimnames attribute from geno_coded.
#dimnames(geno_coded) <- NULL
geno_coded[1:5, 1:5]
# Now check the structure to confirm:
str(geno_coded)

# Compute the additive relationship matrix (GRM) using VanRaden method.
myKI <- sommer::A.mat(geno_coded)  # Compute GRM; min.MAF=0 includes all SNPs.
#myKI <- myKI + diag(1e-4, ncol(myKI), ncol(myKI))  # Add small diagonal value to ensure positive definiteness.

# Optional: Preview and visualize the GRM.
myKI[1:5, 1:5]  # View a subset of the GRM.
colfunc <- colorRampPalette(c("steelblue4", "springgreen", "yellow"))  # Define color palette for heatmap.
png("bw8_Amatrix_updated.png", width = 10, height = 6, units = "in", res = 300)  # Open PNG device for saving plot.
heatmap(myKI, col = colfunc(100), Colv = "Rowv")  # Plot heatmap of GRM (symmetric clustering).
dev.off()  # Close the PNG device.

# Step 5: Prepare genotype Data and phenotype data for GWASTools
# Add index to BIM for SNP IDs and set as row names for genotype matrix.
bim_bw8_corrected <- cbind(index = 1:nrow(bim_bw8_corrected), bim_bw8_corrected)  # Add SNP index column.
all.equal(rownames(X_bw8_subset), bim_bw8_corrected$id)
rownames(X_bw8_subset) <- bim_bw8_corrected$index  # Set row names of genotype matrix to SNP indices.
# Add scanID to phenotype data 
mydat_ordered$scanID <- rownames(mydat_ordered)  # Add a scanID column (row names as IDs; useful for GWASTools).

# Create MatrixGenotypeReader object (handles large genotype data efficiently).
geno <- MatrixGenotypeReader(
  genotype = X_bw8_subset,  # Subsetted genotype matrix.
  snpID = bim_bw8_corrected$index,  # Unique integer IDs for SNPs.
  chromosome = as.integer(bim_bw8_corrected$chr),  # Chromosome codes (as integers; map non-numeric if needed).
  position = as.integer(bim_bw8_corrected$pos),  # Base pair positions.
  scanID = as.integer(mydat_ordered$scanID),  # Unique integer IDs for samples.
  autosomeCode = 1L:33L,  # Custom: Chicken autosomes (1-33).
  XchromCode = 36L,  # Custom: Z chromosome (often X-like in birds).
  YchromCode = 35L,  # Custom: W chromosome (often Y-like in birds).
  XYchromCode = 37L,  # Custom: Pseudo-autosomal if applicable.
  MchromCode = 34L   # Custom: Mitochondrial.
)
genoData <- GenotypeData(geno)  # Wrap in GenotypeData for association tests.

# Create ScanAnnotationDataFrame for phenotypes.
scanAnno <- ScanAnnotationDataFrame(mydat_ordered)  # Convert ordered phenotypes to annotation object.
scanAnno  # Preview the annotation data.

# Step 6: Fit the Null Model
# Identify the phenotype column (here, assuming 'F2_BW8' is the 37th column; verify with match).
columnIndex <- match("F2_BW8", colnames(scanAnno))  # Find column index of the trait.
colnames(scanAnno)[columnIndex]  # Confirm the column name.

# Step 7: Select variables for GWAS model
# Fit a linear mixed model with GRM.
# sire effect?
mix1 <- lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_sire),  # Formula: BW8 as response, random effects for ID (genetic) and dam.
  relmat = list(scanID = myKI),  # Use GRM for genetic relatedness.
  verbose = TRUE,  # Print progress.
  data = mydat_ordered # Use ordered phenotype data.
)

summary(mix1)  # Full model summary.

# Extract variance components and compute heritability.
vc <- VarCorr(mix1)  # Get variance-covariance components.
print(vc, comp = c("Variance"))  # Print variances.
ve <- attr(VarCorr(mix1), "sc")^2  # Residual variance (sigma^2).
h2 <- vc$scanID / (vc$scanID + vc$F1_sire + ve) #the ratio of additive genetic variance (V_A) to total phenotypic variance (V_P)
# not the full narrow-sense h² from all genetic variants (e.g., rare or untyped ones are missed)!
as.numeric(h2)  # Convert to numeric for display.
# sire effect = ~ 4%
vc$F1_sire / (vc$scanID + vc$F1_sire + ve)

# Dam effect?
mix2 <- lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_dam),  # Formula: BW8 as response, random effects for ID (genetic) and dam.
  relmat = list(scanID = myKI),  # Use GRM for genetic relatedness.
  verbose = TRUE,  # Print progress.
  data = mydat_ordered # Use ordered phenotype data.
)

summary(mix2)  # Full model summary.

# Extract variance components and compute heritability.
vc <- VarCorr(mix2)  # Get variance-covariance components.
print(vc, comp = c("Variance"))  # Print variances.
ve <- attr(VarCorr(mix2), "sc")^2  # Residual variance (sigma^2).
h2 <- vc$scanID / (vc$scanID + vc$F1_dam + ve)
as.numeric(h2)  # Convert to numeric for display.
# Dam effect = ~3%
vc$F1_dam / (vc$scanID + vc$F1_dam + ve)


# Sire and Dam effect?
mix3 <- lmebreed(
  F2_BW8 ~ F2_SEX + (1|scanID) + (1|F1_sire) + (1|F1_dam),  # Formula: BW8 as response, random effects for ID (genetic) and dam.
  relmat = list(scanID = myKI),  # Use GRM for genetic relatedness.
  verbose = TRUE,  # Print progress.
  data = mydat_ordered # Use ordered phenotype data.
)

summary(mix3)  # Full model summary.

# Extract variance components and compute heritability.
vc <- VarCorr(mix3)  # Get variance-covariance components.
print(vc, comp = c("Variance"))  # Print variances.
ve <- attr(VarCorr(mix3), "sc")^2  # Residual variance (sigma^2).
h2 <- vc$scanID / (vc$scanID + vc$F1_sire + vc$F1_dam + ve)  # Heritability: genetic var / (genetic + residual var).
as.numeric(h2)  # Convert to numeric for display.
# Dam effect = ~2.8%
vc$F1_dam / (vc$scanID + vc$F1_sire + vc$F1_dam + ve) 
# Sire effect = ~3.5%
vc$F1_sire / (vc$scanID + vc$F1_sire + vc$F1_dam + ve) 

# Set F1_dam as additional random effects 
# Assuming mydat_ordered has the column 'F1_dam' with dam IDs (as character or factor).
str(mydat_ordered)
# Step 1): Create incidence matrix Z for dams (samples x dam levels).
dam_ids <- factor(mydat_ordered$F1_dam)  # Convert to factor to handle levels.
Z <- model.matrix(~ dam_ids - 1)  # Incidence matrix: rows = samples, columns = dam levels, 0/1 entries.
rownames(Z) = mydat_ordered$F2_id
# verified the first column of Z 
# merged5|>filter(F2_id == "B1-11") |>distinct(F1_dam) etc...

# Step 2): Compute dam covariance matrix K_dam = Z %*% t(Z).
# This is efficient; for large n, consider sparseMatrix from Matrix package if memory is an issue.
K_dam <- tcrossprod(Z)  # Results in n x n matrix: 1 if same dam, 0 otherwise.

# Optional: Add small diagonal for stability (similar to GRM).
K_dam <- K_dam + diag(1e-4, nrow(K_dam), ncol(K_dam))

# Step 3: Verify dimensions and alignments.
dim(K_dam)  # Should be [n_samples, n_samples], e.g., [1644, 1644].


# Verify alignments.
all.equal(rownames(myKI), colnames(myKI))  # Ensure GRM is square and symmetric.
all.equal(rownames(K_dam), colnames(K_dam))  # Ensure GRM is square and symmetric.
all.equal(colnames(myKI), colnames(K_dam))
all.equal(rownames(myKI), rownames(K_dam))
all.equal(colnames(myKI), colnames(X_bw8_subset))  # Check if GRM columns match subsetted genotypes.
all.equal(colnames(K_dam), colnames(X_bw8_subset))  # Check if dam matrix columns match subsetted genotypes.
all.equal(mydat_ordered$F2_id, colnames(myKI))
all.equal(mydat_ordered$F2_id, colnames(K_dam))
rownames(myKI) <- mydat_ordered$scanID  # Set row names of GRM to scanIDs.
colnames(myKI) <- mydat_ordered$scanID  # Set column names similarly.
rownames(K_dam) <- mydat_ordered$scanID  # Set row names of dam matrix to scanIDs.
colnames(K_dam) <- mydat_ordered$scanID  # Set column names similarly.

# Step 4): Fit the null model with multiple random effects.
nullmod <- fitNullModel(
  scanAnno,  # Phenotype annotations (unchanged).
  outcome = colnames(scanAnno)[columnIndex],  # Trait name (e.g., "F2_BW8").
  covars = "F2_SEX",  # Only fixed effects; remove "F1_dam" since it's now random.
  cov.mat = list(gen = myKI, dam = K_dam),  # List of covariance matrices for random effects.
  family = "gaussian"  # Assume normal distribution for quantitative trait.
)

# Optional: Preview the fitted model to check variance components.
nullmod$varComp  # Shows estimated variances for each component (gen, dam, residual).

# Fit null model with GRM and covariates (F2_SEX as fixed, F1_dam as fixed).
nullmod2 <- fitNullModel(
  scanAnno,  # Phenotype annotations.
  outcome = colnames(scanAnno)[columnIndex],  # Trait name (e.g., "F2_BW8").
  covars = c("F2_SEX", "F1_dam"),  # Covariates: sex and dam as fixed effects.
  cov.mat = myKI,  # GRM to account for relatedness.
  family = "gaussian"  # Assume normal distribution for quantitative trait.
)

# Step 5): Run GWAS Association Test
# *Key Considerations in GWAS: assumes additivity: Ignores dominance/epistasis!
# Create iterator for efficient block-wise processing of SNPs.
genoIterator <- GenotypeBlockIterator(genoData, snpBlock = 10000)  # Process 10,000 SNPs per block.

# Run single-variant association test 
# 1) using dam as random effects
assoc_K_dam <- assocTestSingle(
  genoIterator,  # Genotype iterator.
  null.model = nullmod,  # Fitted null model.
  BPPARAM = BiocParallel::SerialParam()  # Use serial processing (change to MulticoreParam for parallel).
)

# 2) using dam as fixed effects
assoc_fixed_dam <- assocTestSingle(
  genoIterator,  # Genotype iterator.
  null.model = nullmod2,  # Fitted null model.
  BPPARAM = BiocParallel::SerialParam()  # Use serial processing (change to MulticoreParam for parallel).
)

# Step 8: Adjust P-Values and Filter Significant SNPs
# 1) dam as random effect
assoc_K_dam$p_adj.BH <- p.adjust(assoc_K_dam$Score.pval, method = "BH")  # Benjamini-Hochberg FDR adjustment.
assoc_K_dam$p_adj.bonferroni <- p.adjust(assoc_K_dam$Score.pval, method = "bonferroni")  # Bonferroni correction.
assoc_sig_K_dam <- assoc_K_dam %>% filter(p_adj.BH < 0.05)  # Filter SNPs with FDR < 0.05.
# Join with BIM to add SNP details (e.g., id, chr) on both chromosome and position.
# bim_bw8_corrected |> distinct(pos)

assoc_sig_K_dam <- left_join(assoc_sig_K_dam, bim_bw8_corrected, by = c("variant.id" = "index"))

# 2) dam as fixed effect
assoc_fixed_dam$p_adj.BH <- p.adjust(assoc_fixed_dam$Score.pval, method = "BH")  # Benjamini-Hochberg FDR adjustment.
assoc_fixed_dam$p_adj.bonferroni <- p.adjust(assoc_fixed_dam$Score.pval, method = "bonferroni")  # Bonferroni correction.
assoc_sig_fixed_dam <- assoc_fixed_dam %>% filter(p_adj.BH < 0.05)  # Filter SNPs with FDR < 0.05.
# Join with BIM to add SNP details (e.g., id, chr) on both chromosome and position.
# bim_bw8_corrected |> distinct(pos)

assoc_sig_fixed_dam <- left_join(assoc_sig_fixed_dam, bim_bw8_corrected, by = c("variant.id" = "index"))

# Step 9: Compare with External GCTA Results (Optional)
# Assume 'final_gwas_results_annotated' has GCTA results.
gcta_mlma_sig_SNPs <- final_gwas_results_annotated |> filter(gwas_trait == "bw8") |> distinct(id)  # GCTA sig SNPs.
assoc_sig_K_dam_SNPs <- assoc_sig_K_dam |> distinct(id)  # R GWAS sig SNPs.
common_snps <- intersect(gcta_mlma_sig_SNPs$id, assoc_sig_K_dam_SNPs$id)  # Find overlapping significant SNPs.

assoc_sig_fixed_dam_SNPs <- assoc_sig_fixed_dam |> distinct(id)  # R GWAS sig SNPs.
common_snps2 <- intersect(gcta_mlma_sig_SNPs$id, assoc_sig_fixed_dam_SNPs$id)  # Find overlapping significant SNPs.

# Optional: Save results.
write.csv(assoc_sig_K_dam, "bw8_gwas_results.csv")  # Save significant associations.

# Step 10: QQ plot (with lambda) and circular Manhattan plot
#-----------------------------------------------------------

# install.packages("CMplot")  # run once if CMplot not installed
# library(qqman)
library(CMplot)

## --- QQ plot with genomic inflation factor --------------------------
# Choose the association result you want to visualise:
pvals <- assoc_K_dam$Score.pval      # or assoc_fixed_dam$Score.pval

chisq  <- qchisq(1 - pvals, df = 1)                      # χ² statistics
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)   # inflation factor

qq(pvals,
   main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),
   xlab = "Expected -log10(p)",
   ylab = "Observed -log10(p)")

# Save QQ plot if desired:
png("bw8_gwas_qqplot.png", width = 6, height = 6, units = "in", res = 300)
qq(pvals, main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"))
dev.off()

# OR use CMplot function
# QQ plot (using CMplot) with title including lambda
manplot <- assoc_K_dam |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)

names(manplot)

CMplot(
  Pmap = manplot,
  plot.type = "q",  # QQ plot
  conf.int = TRUE,  # Show confidence intervals
  box = TRUE,  # Add box around plot
  main = paste0("QQ Plot (λ = ", round(lambda, 3), ")"),  # Add title with lambda
  file.output = TRUE,
  file = "tiff"
)


## --- Circular Manhattan plot ---------------------------------------
# Prepare data frame: SNP ID, chromosome, position, p-value
manplot <- assoc_K_dam |>
  dplyr::select(SNP = variant.id, Chromosome = chr, Position = pos, P = Score.pval)

names(manplot)
# Extract raw p-values
raw_p <- manplot$P

# Compute BH-adjusted p-values
adj_p <- p.adjust(raw_p, method = "BH")

# Find the BH threshold: max raw p where adjusted p <= 0.05 (FDR level)
# If none meet the criterion, threshold will be 0 (no significant SNPs)
bh_threshold <- ifelse(any(adj_p <= 0.05), max(raw_p[adj_p <= 0.05], na.rm = TRUE), 0)

# Print the threshold for verification
print(bh_threshold)

# Draw circular Manhattan plot (Bonferroni line shown)
CMplot(
  manplot,
  plot.type   = "c",                   # circular
  LOG10       = TRUE,
  threshold = c(bh_threshold, 0.05 / nrow(manplot)),
  threshold.col = c("orange","red"),
  chr.den.col = NULL,
  file        = "jpg",
  dpi         = 300,
  width       = 9,
  height      = 9
)

# Use the script: manhattan_plots.R (https://github.com/seanjosephjurgens/UKBB_200KWES_CVD/blob/v1.2/manhattan_plots.R)
source("manhattan_plots.R")
library(calibrate)
library(qqman)
CMplot(
  Pmap = manplot,  # Your data frame
  plot.type = "c",  # "c" for circular Manhattan; use "m" for rectangular, "q" for QQ, "d" for density
  col = c("blue", "red"),  # Colors for chromosomes
  threshold = 1e-5,  # Significance threshold line (e.g., Bonferroni)
  threshold.col = "red",  # Color for threshold line
  amplify = TRUE,  # Amplify points above threshold
  signal.cex = 1.2,  # Size of significant points
  file.output = TRUE,  # Save to file (e.g., jpg/pdf)
  file = "jpg",  # Output format
  dpi = 300,  # Resolution
  verbose = TRUE  # Print progress
)

# QQ plot (using CMplot)
CMplot(
  Pmap = manplot,
  plot.type = "q",  # QQ plot
  conf.int = TRUE,  # Show confidence intervals
  box = TRUE,  # Add box around plot
  file.output = TRUE,
  file = "pdf"
)

# SNP Density Plot (using CMplot)
CMplot(
Pmap =  manplot,
plot.type = "d",  # Density plot
bin.size = 1e6,  # Window size for density
file.output = TRUE,
file = "jpg"
)

# All Plot Types at Once
CMplot(
  Pmap = manplot,
  amplify = TRUE,
  threshold = c(bh_threshold, 0.05 / nrow(manplot)),
  threshold.col = c("red","darkred"),
  plot.type=c("m","q","c","d"), # "b" = all (circular + rectangular + QQ + density)
)


### Optional: for multi-traits 
# Assume Pmap has multiple traits (e.g., columns P_Trait1, P_Trait2)
# Convert Chromosome to character for safe replacement
manplot$Chromosome <- as.character(manplot$Chromosome)

# Map specific non-numeric labels to numbers
manplot$Chromosome[manplot$Chromosome == "M"] <- "34"
manplot$Chromosome[manplot$Chromosome == "X"] <- "36"
# Add more if needed, e.g.:
# manplot$Chromosome[manplot$Chromosome == "Y"] <- "35"
# manplot$Chromosome[manplot$Chromosome == "MT"] <- "34"  # If "MT" exists

# Convert to numeric
manplot$Chromosome <- as.numeric(manplot$Chromosome)

# Verify: should all be numeric now
str(manplot$Chromosome)
#manplot$pch <- 19
# Rerun the plot function
multi_manhattan(
  x = manplot, 
  chr = "Chromosome", bp = "Position", p = "P", snp = "SNP",
  col = c("gray10", "gray60"),
  suggestiveline = -log10(1e-5),
  genomewideline = -log10(5e-8),
  logp = TRUE,
  file.output = TRUE,
  file = "png"
)

## Add gene annotation for the significant SNPs!
#Step 1: Install and Load Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))
install.packages("dplyr")  # For data manipulation

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

#Step 2: Import the GTF File
# Path to your GTF file
gtf_file <- "Gallus_gallus_gca000002315v5.GRCg6a.114.gtf.gz"

# Import GTF
gtf <- import(gtf_file)

# Filter to genes only (or exons/transcripts if preferred; adjust as needed)
genes <- gtf[gtf$type == "gene"]  # GRanges object with gene metadata

# Preview
head(genes)
#Step 3: Prepare Your Significant SNPs
# Example: Your significant SNPs data frame (replace with yours)
# assoc_sig <- data.frame(
#   id = c("SNP1", "SNP2"),
#   chr = c(1, 2),  # Chromosome (match GTF format, e.g., "1" or "chr1"; check unique(genes$seqnames))
#   pos = c(1000000, 2000000)
# )

# Define window size (bp; set to 0 for exact position overlap only) 
window_size <- 1000000

# Replace "X" with "Z" in seqnames (convert to character first if needed)
assoc_sig_K_dam <- assoc_sig_K_dam %>%
  mutate(chr.x = as.character(chr.x)) %>%  # Ensure it's character
  mutate(chr.x = ifelse(chr.x == "X", "Z", chr.x))
# Replace "M" with "MT" in seqnames
assoc_sig_K_dam <- assoc_sig_K_dam %>%
  mutate(chr.x = as.character(chr.x)) %>%  # Ensure it's character
  mutate(chr.x = ifelse(chr.x == "M", "MT", chr.x))

# Verify the change
assoc_sig_K_dam %>% distinct(chr.x)


# Create GRanges for SNPs with optional windows
snp_ranges <- makeGRangesFromDataFrame(
  assoc_sig_K_dam %>%
    mutate(
      start = pos.x - window_size,
      end = pos.x + window_size,
      strand = "*"  # Assume unstranded
    ),
  keep.extra.columns = TRUE,  # Retain SNP metadata (id, etc.)
  na.rm = TRUE,
  seqnames.field = "chr.x"  # Use your chromosome column
)

# Step 4: Find Overlapping Genes
# Use findOverlaps to identify genes overlapping each SNP (or window).
# Find overlaps between SNPs and genes
overlaps <- findOverlaps(snp_ranges, genes, type = "any")  # "any" for partial overlap; use "within" for SNPs inside genes

# Summarize genes per SNP
gene_summary <- data.frame(
  snp_index = queryHits(overlaps),
  gene_index = subjectHits(overlaps)
) %>%
  mutate(
    snp_id = assoc_sig_K_dam$id[snp_index],
    gene_id = mcols(genes)$gene_id[gene_index],
    gene_name = mcols(genes)$gene_name[gene_index],
    gene_biotype = mcols(genes)$gene_biotype[gene_index],
    gene_start = start(genes)[gene_index],
    gene_end = end(genes)[gene_index],
    gene_strand = as.character(strand(genes)[gene_index])
  ) %>%
  group_by(snp_id) %>%
  summarise(
    genes = paste(unique(gene_name), collapse = "; "),
    details = paste(gene_id, " (", gene_start, "-", gene_end, ", strand: ", gene_strand, ", biotype: ", gene_biotype, ")", collapse = "; ")
  )

# Join to original SNPs and handle no-overlap cases
annotated_snps <- assoc_sig_K_dam %>%
  left_join(gene_summary, by = c("id" = "snp_id")) %>%
  mutate(genes = ifelse(is.na(genes), "No genes in window", genes))

# Preview results
head(annotated_snps)

# Export to CSV
write.csv(annotated_snps, "annotated_significant_snps_bw8_window1Mb.csv", row.names = FALSE)