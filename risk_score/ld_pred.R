
##
#
# This script runs LDpred2 on the MoBa genotypes.
#
# Notes: 
# - It is better to avoid running this script in Workbench due to the high memory requirements.
# - In order to adapt this script to your summary stats, edit the path to 'score file', the 'score_name', and GWAS summary stats parsing.
# - The script stores large intermediate files in the output folder, you can delete these once you are done with the score computation.
# - The script uses a lot of memory, be mindful of other users when computing multiple scores in parallel.
#
##

# Analysis files, need to be updated for every analysis

args <- commandArgs(TRUE)

score_file <- args[1]
score_name <- args[2]

documentation_folder <- args[3] # For QC and documentation, no individual-level data should be saved there

output_folder <- args[4]
intermediate_folder <- file.path(output_folder, "ldpred2_files")


# Housekeeping

if (!file.exists(score_file)) {
  
  stop(paste0("GWAS summary stats file ", score_file, " not found."))
  
}
if (!dir.exists(documentation_folder)) {
  
  dir.create(documentation_folder, recursive = T)
  
}
if (!dir.exists(intermediate_folder)) {
  
  dir.create(intermediate_folder, recursive = T)
  
}


# MoBa Files, do not need to be changed

bed_file <- "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2025.09.25/moba_genotypes_2025.09.25_common.bed"
backing_file_base <- "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2025.09.25/ldpred2/moba_genotypes_2025.09.25_common"


# Libraries

library(bigsnpr)
library(janitor)
library(glue)
library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(grid)


# Parameters

set.seed(20260210)
n_cores <- 8

options(bigstatsr.check.parallel.blas = FALSE)


# Set up genotypes parser

print(glue("{Sys.time()} - Setting up genotypes parser"))

obj.bigSNP <- snp_attach(glue("{backing_file_base}.rds"))

G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))


# LD reference

print(glue("{Sys.time()} - Reading LD reference map"))

map_ldref <- readRDS("/mnt/archive/hapmap/map_hm3_ldpred2.rds")


# Read GWAS summary statistics used for score training
# Note that columns need to contain 'chr', 'pos', 'a0', 'a1', 'beta', 'beta_se', 'n_eff'.

print(glue("{Sys.time()} - Reading summary statistics"))

df_beta <- read.table(
  file = score_file,
  header = T,
  sep = "\t"
) %>% 
  clean_names() %>% 
  select(
    rsid = marker_name,
    chr = chr,
    pos = pos,
    a0 = allele1,
    a1 = allele2,
    beta = effect,
    beta_se = std_err,
    p_value = p_value,
    n_eff = totalsamplesize
  ) %>% 
  mutate(
    a0 = toupper(a0),
    a1 = toupper(a1)
  )

n_variants <- nrow(df_beta)

print(glue("{Sys.time()} - {n_variants} variants loaded"))


# Remove variants with low sample size

n_min <- 0.7*max(df_beta$n_eff)

df_beta <- df_beta %>% 
  filter(
    n_eff > n_min
  )

if (nrow(df_beta) < 10000) {
  
  stop("Less than 10000 variants available after filtering for sample size, please verify why the variants get filtered out.")
  
}

print(glue("{Sys.time()} - {n_variants - nrow(df_beta)} removed due to low sample size, {nrow(df_beta)} variants remaining"))
n_variants <- nrow(df_beta)


# Match with MoBa

print(glue("{Sys.time()} - Matching with MoBa"))

df_beta <- snp_match(df_beta, map)

if (nrow(df_beta) < 10000) {
  
  stop("Less than 10000 variants available after matching with MoBa, please verify why the variants get filtered out.")
  
}

print(glue("{Sys.time()} - {n_variants - nrow(df_beta)} removed after matching with MoBa, {nrow(df_beta)} variants remaining"))
n_variants <- nrow(df_beta)


# Match with ldref

print(glue("{Sys.time()} - Matching with LD ref"))

df_beta <- df_beta %>% filter(
    chr != "X" & chr != "Y"
  ) %>%
  mutate(
    chr = as.numeric(chr)
  )

df_beta <- snp_match(df_beta, map_ldref)

if (nrow(df_beta) < 10000) {
  
  stop("Less than 10000 variants available after matching with LD ref, please verify why the variants get filtered out.")
  
}

print(glue("{Sys.time()} - {n_variants - nrow(df_beta)} removed after matching to LD ref, {nrow(df_beta)} variants remaining"))
n_variants <- nrow(df_beta)


# Quality control of summary statistics as described here: https://www.sciencedirect.com/science/article/pii/S2666247722000525

print(glue("{Sys.time()} - Quality control of the summary stats"))

sd_ldref <- with(df_beta, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_y <- with(df_beta, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01)))
sd_ss <-  with(df_beta, sd_y / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

sd_mismatch <- qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

png(
  filename = glue("{documentation_folder}/snp_qc_plot.png"),
  width = 800,
  height = 600
)
grid.draw(sd_mismatch)
dev <- dev.off()

df_beta <- df_beta[!is_bad, ]

if (nrow(df_beta) < 10000) {
  
  stop("Less than 10000 variants available after QC, please verify why the variants get filtered out.")
  
}

print(glue("{Sys.time()} - {n_variants - nrow(df_beta)} removed after QC, {nrow(df_beta)} variants remaining"))
n_variants <- nrow(df_beta)


# Correlation matrix

matrix_file <- file.path(intermediate_folder, glue("hm3_corr_{score_name}.sbk"))
corr_file <- file.path(intermediate_folder, glue("hm3_corr_{score_name}.rds"))

if (!file.exists(corr_file)) {
  
  print("{Sys.time()} - Creating sparse SNP correlation matrix")
  
  for (chr in 1:22) {
    
    print(glue("chr {chr} of 22"))
    
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    
    corr_chr <- readRDS(paste0("/mnt/scratch/ldpred2/ld_reference/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
    
    if (chr == 1) {
      corr <- as_SFBM(corr_chr, backingfile = file.path(intermediate_folder, glue("hm3_corr_{score_name}")), compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  saveRDS(corr, corr_file)
  
} else {
  
  corr <- readRDS(corr_file)
  
}


# Heritability estimation from LDSC

print(glue("{Sys.time()} - Estimating heritability"))

ldsc <- with(
  df_beta, 
  snp_ldsc(
    ld, 
    ld_size = nrow(map_ldref),
    chi2 = (beta / beta_se)^2,
    sample_size = n_eff,
    ncores = n_cores
  )
)

h2_est <- ldsc[["h2"]]

print(glue("Heritablity-estimate: {round(h2_est, digits = 3)}"))


# LDpred2-auto

print(glue("{Sys.time()} - Running ldpred2 auto"))

shrinkage_coeff <- 0.95

rds_file <- file.path(intermediate_folder, "multi_auto.rds")

if (!file.exists(rds_file)) {
  
  multi_auto <- snp_ldpred2_auto(
    corr = corr, 
    df_beta = df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.9, 30),
    ncores = n_cores, 
    shrink_corr = shrinkage_coeff,
    allow_jump_sign = TRUE,
    use_MLE = TRUE
  )
  
  saveRDS(multi_auto, file = rds_file)
  
} else {
  
  multi_auto <- readRDS(rds_file)
  
}

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

if (sum(!is.na(beta_auto)) == 0) {
  
  print("{Sys.time()} - All 'beta_est' values are NA after 'snp_ldpred2_auto', rerunning with 'use_MLE' turned off.")
  
  multi_auto <- snp_ldpred2_auto(
    corr = corr, 
    df_beta = df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.9, 30),
    ncores = n_cores, 
    shrink_corr = shrinkage_coeff,
    allow_jump_sign = TRUE,
    use_MLE = FALSE
  )
  
  saveRDS(multi_auto, file = rds_file)
  
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  
}


# Filter and combine based on correlation estimates

print(glue("{Sys.time()} - Filtering runs and combine based on correlation estimates"))

range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
beta_auto_avg <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))


# Check chains for convergence, create common plot

print(glue("{Sys.time()} - Checking chains for convergence and making QC plots"))

auto_params <- rbindlist(lapply(multi_auto, function(x) {
  data.table(
    p_init = x$p_init, 
             h2_init = x$h2_init, 
             p_est = x$p_est, 
             h2_est = x$h2_est
    )
}))

auto_params <- auto_params %>%
  mutate(
    paramset = as.numeric(rownames(.)),
    keep = ifelse(paramset %in% keep, "chain kept", "chain excluded")
  )

auto_path <- foreach(pIdx = seq_along(multi_auto), .combine=rbind) %do% {
  auto = multi_auto[[pIdx]]
  data.table(paramset = pIdx, path_iter = seq_along(auto$path_p_est), 
             p_est = auto$path_p_est, h2_est = auto$path_h2_est)
}

p_plot <- ggplot(auto_path) + aes(x = path_iter, y=p_est) +
  theme_bigstatsr() + 
  geom_hline(data = auto_params, aes(yintercept=p_est, colour = keep)) +
  geom_point(shape=19, size=0.5) +
  scale_y_log10(name="p") + xlab("") +
  facet_wrap(~ paramset, ncol=10, labeller = label_both) + 
  theme(strip.background=element_blank(), strip.text=element_text(size=6), 
        axis.text=element_text(size=6), axis.title=element_text(size=10),
        legend.position = "top", 
        legend.title = element_blank())

png(
  filename = file.path(documentation_folder, "p_est_qc.png"),
  width = 900,
  height = 900
)

grid.draw(p_plot)

dev <- dev.off()


h2_plot <- ggplot(auto_path) + aes(x = path_iter, y=h2_est) +
  theme_bigstatsr() + 
  geom_hline(data = auto_params, aes(yintercept=h2_est, colour = keep)) +
  geom_point(shape=19, size=0.5) +
  ylab("h2") + xlab("") +
  facet_wrap(~ paramset, ncol=10, labeller = label_both) +
  theme(strip.background=element_blank(), 
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=10),
        legend.position = "none")

png(
  filename = file.path(documentation_folder, "h2_est_qc.png"),
  width = 900,
  height = 900
)

grid.draw(h2_plot)

dev <- dev.off()


# Write weights to plink-friendly format

print(glue("{Sys.time()} - Exporting results"))

weights_plink <- data.frame(
  rsid = df_beta$rsid,
  allele1 = df_beta$a1,
  weight = beta_auto_avg
)

write.table(
  weights_plink, 
  file = file.path(output_folder, glue("{score_name}_weights_plink_format")),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)

  

