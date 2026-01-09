
##
#
# This script prepares phenotype files for Regenie
#
##

# Seed for random choice of sample

set.seed(12122025)

# Libraries

library(conflicted)
library(janitor)
library(glue)
library(dplyr)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Command line input
  
args <- commandArgs(TRUE)

if (length(args) != 4) {
  
  stop(paste0("Four command line arguments expected. ", length(args), " found."))
  
}

pheno_table <- args[1]

if (!file.exists(pheno_table)) {
  
  stop(paste0("Phenotype table ", pheno_table, " not found."))
  
}

psam_file <- args[2]

if (!file.exists(psam_file)) {
  
  stop(paste0("Psam file ", psam_file, " not found."))
  
}

pcs_file <- args[3]

if (!file.exists(pcs_file)) {
  
  stop(paste0("PCs file ", pcs_file, " not found."))
  
}

gwas_pheno_folder <- args[4]

if (!dir.exists(gwas_pheno_folder)) {
  
  stop(paste0("GWAS phenotypes folder ", gwas_pheno_folder, " not found."))
  
}


# Pheno file

print(paste0(Sys.time(), " - Loading phenotypes from '", pheno_table, "'"))

pheno_table <- read.table(
  file = pheno_table,
  header = T,
  sep = "\t"
)

required_columns <- c("child_sentrix_id", "mother_sentrix_id", "father_sentrix_id", "child_batch", "mother_batch", "father_batch")

for (column in required_columns) {
  
  if (!column %in% names(pheno_table)) {
    
    stop(paste0("Column '", column, "' not found in `", pheno_table, "`."))
    
  }
}


# Psam file

print(paste0(Sys.time(), " - Loading psam file"))

psam_df <- read.table(
  file = psam_file,
  header = F,
  sep = "\t",
  stringsAsFactors = F
)
psam_df <- psam_df[, 1:2]
names(psam_df) <- c("FID", "IID")


# Pcs file

print(paste0(Sys.time(), " - Loading pcs file"))

pcs_df <- read.table(
  file = pcs_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  select(
    "IID", starts_with("PC")
  )



# Set up data frames for GWAS

print(paste0(Sys.time(), " - Setting up for gwas"))

pheno_table_gwas_child <- psam_df %>% 
  inner_join(
    pheno_table %>%
      filter(
        !is.na(child_sentrix_id)
      ) %>% 
      select(
        IID = child_sentrix_id,
        child_batch,
        where(is.numeric)
      ),
    by = "IID",
    multiple = "any"
  ) %>% 
  left_join(
    pcs_df,
    by = "IID"
  )

for (batch_name in unique(pheno_table_gwas_child$child_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_child[[batch_column]] <- ifelse(pheno_table_gwas_child$child_batch == batch_column, 1, 0)
  
}

pheno_table_gwas_mother <- psam_df %>% 
  inner_join(
    pheno_table %>% 
      filter(
        !is.na(mother_sentrix_id)
      ) %>% 
      select(
        IID = mother_sentrix_id,
        mother_batch,
        where(is.numeric)
      ),
    by = "IID",
    multiple = "any"
  ) %>% 
  left_join(
    pcs_df,
    by = "IID"
  )

for (batch_name in unique(pheno_table_gwas_mother$mother_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_mother[[batch_column]] <- ifelse(pheno_table_gwas_mother$mother_batch == batch_column, 1, 0)
  
}

pheno_table_gwas_father <- psam_df %>% 
  inner_join(
    pheno_table %>% 
      filter(
        !is.na(father_sentrix_id)
      ) %>% 
      select(
        IID = father_sentrix_id,
        father_batch,
        where(is.numeric)
      ),
    by = "IID",
    multiple = "any"
  ) %>% 
  left_join(
    pcs_df,
    by = "IID"
  )

for (batch_name in unique(pheno_table_gwas_father$father_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_father[[batch_column]] <- ifelse(pheno_table_gwas_father$father_batch == batch_column, 1, 0)
  
}

pheno_table_gwas_parents_mother <- pheno_table_gwas_mother

for (column in names(pheno_table_gwas_father)[!pheno_table_gwas_father %in% pheno_table_gwas_mother]) {
  
  pheno_table_gwas_parents_mother[[column]] <- 0
  
}

pheno_table_gwas_parents_father <- pheno_table_gwas_father

for (column in names(pheno_table_gwas_mother)[!pheno_table_gwas_mother %in% pheno_table_gwas_father]) {
  
  pheno_table_gwas_parents_father[[column]] <- 0
  
}

pheno_table_gwas_parents <- rbind(pheno_table_gwas_parents_mother, pheno_table_gwas_parents_father)


# Write tables

print(paste0(Sys.time(), " - Export of tables to ", gwas_pheno_folder))

write.table(
  x = pheno_table_gwas_child,
  file = file.path(gwas_pheno_folder, "pheno_children"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_mother,
  file = file.path(gwas_pheno_folder, "pheno_mothers"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_father,
  file = file.path(gwas_pheno_folder, "pheno_fathers"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_father,
  file = file.path(gwas_pheno_folder, "pheno_parents"),
  row.names = F,
  col.names = T,
  quote = F
)
