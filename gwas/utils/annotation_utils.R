##
#
# This script contains functions for the query of Phenoscanner and Ensembl. 
# 
##


# Ensembl genes

genes <- read.table(
  "gwas/utils/gene_coordinates.GRCh37.87.gz",
  header = T,
  sep = "\t"
)


#' Fixes custom ids to make valid file names.
#' 
#' @param variantId the id of the variant
#' 
#' @return the cleaned id of the variant
cleanRsId <- function(
    variantIds
) {
  
  # check for rs ids with a suffix
  
  suffixI <- !is.na(variantIds) & is.character(variantIds) & startsWith(variantIds, "rs") & str_detect(variantIds, '_')
  index <- ifelse(!is.na(variantIds) & is.character(variantIds), str_locate(variantIds, pattern = "_")[, 1], 0)
  
  cleanIds <- ifelse(!is.na(variantIds), variantIds, "")
  cleanIds <- ifelse(is.character(variantIds), variantIds, as.character(variantIds))
  cleanIds <- ifelse(suffixI, str_sub(variantIds, end = index - 1), variantIds)
  
  # Remove special characters 
  
  cleanIds <- str_replace_all(cleanIds, "/", "_")
  cleanIds <- str_replace_all(cleanIds, ":", "_")
  
  return(cleanIds)
  
}


#' Writes documentation for the nearest genes.
#' 
#' @param variantId the id of the variant to query
#' @param chromosome the chromosome of the variant to query
#' @param bp the pb of the variant to query
#' @param maxDistance the distance to use around variants of interest for gene annotation
#' 
#' @return returns a label for the table or an empty label
get_nearest_gene_docs <- function(
    variantId,
    chromosome,
    bp,
    ensemblFolder,
    maxDistance = 250000
) {
  
  nearest_genes <- genes %>% 
    filter(
      chr == chromosome & bp <= end + maxDistance & bp >= start - maxDistance
    )
  
  if (nrow(nearest_genes) == 0) {
    
    return("No gene found")
    
  }
  
  cleanVariantId <- cleanRsId(variantId)
  
  ensemblFile <- file.path(ensemblFolder, paste0(cleanVariantId, ".md"))
  
  write(x = glue("# {variantId}\n"), file = ensemblFile, append = F)
  
  nearest_genes <- nearest_genes %>%
    distinct() %>% 
    mutate(
      distance = case_when(
        bp < start ~ start - bp,
        bp > end ~ bp - end,
        T ~ 0
      ),
      nearest_gene = ifelse(distance == min(distance), T, F)
    )
  
  write(x = "## Closest genes\n", file = ensemblFile, append = T)
  
  write(x = glue("Genes within {maxDistance} bp from SNP (chr {chromosome}, pos {bp})\n"), file = ensemblFile, append = T)
  
  
  write(x = paste0("| Chr | Start | End | Name | Distance |"), file = ensemblFile, append = T)
  write(x = paste0("| --- | ----- | --- | ---- | -------- |"), file = ensemblFile, append = T)
  
  for (i in 1:nrow(nearest_genes)) {
    
    if (nearest_genes$nearest_gene[i]) {
      
      gene_prefix <- "**_"
      gene_suffix <- "**_"
      
    } else {
      
      gene_prefix <- "_"
      gene_suffix <- "_"
      
    }
    
    write(x = paste0("| ", nearest_genes$chr[i], " | ", nearest_genes$start[i], " | ", nearest_genes$end[i], " | ", gene_prefix, nearest_genes$gene_name[i], gene_suffix, " | ", nearest_genes$distance[i], " |"), file = ensemblFile, append = T)
    
  }
  
  write(x = "\n", file = ensemblFile, append = T)

  nearest_gene_name <- paste(nearest_genes$gene_name[nearest_genes$nearest_gene], collapse = ", ")
  
  return(nearest_gene_name)
  
}
