
##
#
# This merges, formats, and documents regenie GWAS results.
#
##



#' Returns a data frame with the to hits.
#' 
#' @param association_df the association data frame with results from regenie
#' 
#' @return a data frame containing the association results for the top hits
get_top_hits <- function(
    association_df
) {

  temp_df <- association_df %>% 
    filter(
      log10p > -log10(p_value_threshold)
    )
  
  annotation_df <- list()
  
  while (nrow(temp_df) > 0) {
    
    top_snp_i <- which.max(temp_df$log10p)
    
    annotation_df[[length(annotation_df) + 1]] <- temp_df[top_snp_i, ]
    
    bp_min <- temp_df$genpos[top_snp_i] - bp_limit
    bp_max <- temp_df$genpos[top_snp_i] + bp_limit
    
    temp_df <- temp_df %>% 
      filter(
        genpos < bp_min | genpos > bp_max
      )
    
  }
  
  if (length(annotation_df) > 0) {
    
    annotation_df <- do.call(rbind, annotation_df)
    
  } else {
    
    annotation_df <- temp_df
    
  }
  
  return(annotation_df)
  
}


#' Returns a ggplot object with the MH for the given association data frame.
#' 
#' @param association_df the association data frame
#' @param association_df the data frame with variants to annotate
#' 
#' @return a ggplot object with the MH
get_mh <- function(
    association_df,
    annotation_df
) {
  
  # Arrange for plotting
  
  plot_df <- association_df %>%
    arrange(
      log10p
    ) %>%
    mutate(
      chromosome_number = as.numeric(ifelse(chrom == 'X', 23, chrom)),
      x = chromosome_start[chromosome_number] + genpos,
      color = chromosome_number %% 2
    ) %>%
    arrange(
      chromosome_number, log10p, genpos
    )
  
  maxP <- max(plot_df$log10p)
  
  colors <- c(mh_color_1, mh_color_2)
  
  if (!is.null(annotation_df) && nrow(annotation_df) > 0) {
    
    plot_df$label <- NA
    
    plotannotation_df <- annotation_df %>% 
      mutate(
        chromosome_number = as.numeric(ifelse(chrom == 'X', 23, chrom)),
        x = chromosome_start[chromosome_number] + genpos
      ) %>% 
      arrange(
        x
      )
    
    for (i in 1:nrow(plotannotation_df)) {
      
      id <- plotannotation_df$id[i]
      chrom <- plotannotation_df$chrom[i]
      genpos <- plotannotation_df$genpos[i]
      
      if (is.na(id)) {
        
        stop("NA snp")
        
      }
      
      plot_df$label[plot_df$id == id] <- id
      plot_df$color[plot_df$id == id] <- 1 + i
      
      plot_df$color[plot_df$chrom == chrom & plot_df$genpos >= genpos - bp_limit & plot_df$genpos <= genpos + bp_limit] <- 1 + i
      
    }
    
    colors <- c(
      colors, 
      scico(
        n = nrow(annotation_df), 
        palette = "batlow", 
        end = 0.8
      )
    )
  }
  
  plot_df <- plot_df %>% 
    mutate(
      color = factor(color, levels = 0:(length(colors)-1))
    ) %>% 
    arrange(
      color
    )
  
  # Chromosome labels
  
  x_labels <- 1:22
  x_labels[x_labels %% 2 == 0 & x_labels > 17] <- ""
  x_labels <- c(x_labels, "X")
  
  # y axis
  
  max_y <- 5 * ceiling(max(plot_df$log10p / 5))
  max_y <- max(max_y, 10)
  
  y_breaks <- c(0, 5, -log10(5e-8))
  y_labels <- c("", 5, round(-log10(5e-8), digits = 1))
  
  last_break <- floor(max(max_y / 10))
  
  if (last_break > 0) {
    
    newBreaks <- 10*(1:last_break)
    
    while(length(newBreaks) > 3) {
      
      newBreaks <- newBreaks[c(T, F)]
      
    }
    
    y_breaks <- c(y_breaks, newBreaks)
    y_labels <- c(y_labels, round(newBreaks, digits = 1))
    
  }
  
  
  # Build plot
  mh_plot <- ggplot() + 
    geom_hline(
      yintercept = -log10(5e-8), 
      col = "green4", 
      linewidth = 0.3
    )
  
  if (!is.null(annotation_df) && nrow(annotation_df) > 0 && length(unique(annotation_df$snp)) < 20) {
    
    max_y <- max_y + 2
    
    annotatedDF <- plot_df %>% 
      filter(
        !is.na(label)
      )
    
    mh_plot <- mh_plot +
      geom_segment(
        data = annotatedDF,
        mapping = aes(
          x = x,
          xend = x,
          y = log10p,
          yend = 1.05 * maxP,
          col = color
        ),
        linetype = "dotted"
      ) +
      geom_text_repel(
        data = annotatedDF,
        mapping = aes(
          x = x,
          y = 1.05 * maxP,
          label = label,
          col = color
        ),
        direction = "x",
        hjust = 0,
        angle = 90,
        nudge_y = 1
      )
    
  }
  
  mh_plot <- mh_plot + 
    geom_point(
      data = plot_df,
      aes(
        x = x, 
        y = log10p, 
        col = color
      ), 
      size = 2
    ) +
    scale_y_continuous(
      name = "p-value [-log10]", 
      breaks = y_breaks, 
      labels = y_labels, 
      expand = expansion(
        mult = c(0, 0.05)
      ), 
      limits = c(0, max_y)
    ) + 
    scale_x_continuous(
      name = "Chromosome", 
      breaks = chromosome_middle, 
      labels = x_labels, 
      limits = c(0, genome_length), 
      expand = expansion(
        mult = 0.01
      )
    ) + 
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) + 
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3),
      strip.background = element_rect(
        fill = "grey99"
      )
    )
  
  return(mh_plot)
  
}

#' Returns a ggplot object with the QQ for the given association data frame.
#' 
#' @param association_df the association data frame
#' 
#' @return a ggplot object with the plot
get_qq <- function(
    association_df
) {
  
  plot_df <- association_df %>%
    mutate(
      maf_bin = case_when(
        a1freq < 0.001 ~ "< 0.1 %",
        a1freq < 0.01 ~ "0.1-1 %",
        a1freq < 0.1 ~ "1-10 %",
        a1freq >= 0.1 ~ "> 10 %"
      )
    ) %>% 
    group_by(
      maf_bin
    ) %>% 
    arrange(
      log10p
    ) %>% 
    mutate(
      observed_log_p = log10p,
      expected_log_p = sort(-log10(ppoints(n = n())))
    ) %>% 
    ungroup()
  
  # Color
  
  plot_df <- plot_df %>% 
    mutate(
      maf_bin = factor(maf_bin, levels = c("< 0.1 %", "0.1-1 %", "1-10 %", "> 10 %"))
    ) %>% 
    arrange(
      maf_bin
    )
  
  # Axis
  
  max_value <- max(plot_df$observed_log_p, plot_df$expected_log_p, 10)
  
  
  # Make plot
  qq_plot <- ggplot(
    data = plot_df
  ) + 
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted"
    ) + 
    geom_hline(
      yintercept = -log10(5e-8), 
      col = "green4", 
      linewidth = 0.3
    )  +
    geom_point(
      mapping = aes(
        x = expected_log_p,
        y = observed_log_p,
        col = maf_bin
      ),
      size = 2
    ) +
    scale_color_manual(
      values = scico(
        n = 4,
        begin = 0.2,
        end = 0.8
      )
    ) +
    scale_x_continuous(
      name = "Expected p-value [-log10]",
      limits = c(0, max_value),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "p-value [-log10]", 
      limits = c(0, max_value),
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(qq_plot)
  
}

#' Returns a ggplot object with the beta plotted against the maf.
#' 
#' @param association_df the association data frame
#' @param annotation_df the annotation data frame
#' 
#' @return a ggplot object with the plot
get_beta_maf <- function(
    association_df,
    annotation_df
) {
  
  plot_df <- association_df %>%
    filter(
      is.finite(beta) & is.finite(se) & is.finite(a1freq)
    ) %>% 
    arrange(
      a1freq
    )
  
  annotated_points_df <- plot_df %>% 
    filter(
      id %in% annotation_df$id
    )
  
  
  # Make plot
  beta_maf_plot <- ggplot() + 
    geom_hline(
      yintercept = 0
    ) +
    geom_segment(
      data = plot_df,
      mapping = aes(
        x = 100 * a1freq,
        xend = 100 * a1freq,
        y = beta - qnorm(0.975) * se,
        yend = beta + qnorm(0.975) * se
      ),
      col = "grey90",
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_df,
      mapping = aes(
        x = 100 * a1freq,
        y = beta
      ),
      col = "grey70",
      size = 0.8
    ) +
    geom_segment(
      data = annotated_points_df,
      mapping = aes(
        x = 100 * a1freq,
        xend = 100 * a1freq,
        y = beta - qnorm(0.975) * se,
        yend = beta + qnorm(0.975) * se
      ),
      col = "grey30",
      linewidth = 0.8
    ) +
    geom_point(
      data = annotated_points_df,
      mapping = aes(
        x = 100 * a1freq,
        y = beta
      ),
      col = "grey20",
      size = 0.8
    ) +
    scale_x_continuous(
      name = "Minor Allele Frequency [%]",
      limits = c(0, 50),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "Effect Size Estimate [95% CI]",
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(beta_maf_plot)
  
}

#' Returns a ggplot object with the se plotted against the maf.
#' 
#' @param association_df the association data frame
#' @param annotation_df the annotation data frame
#' 
#' @return a ggplot object with the plot
get_sd_maf_plot <- function(
    association_df,
    annotation_df
) {
  
  plot_df <- association_df %>%
    filter(
      is.finite(beta) & is.finite(se) & is.finite(a1freq)
    ) %>% 
    arrange(
      a1freq
    )
  
  annotated_points_df <- plot_df %>% 
    filter(
      id %in% annotation_df$id
    )
  
  
  # Make plot
  se_maf_plot <- ggplot() + 
    geom_hline(
      yintercept = 0
    ) +
    geom_point(
      data = plot_df,
      mapping = aes(
        x = 100 * a1freq,
        y = se
      ),
      col = "grey70",
      size = 0.8
    ) +
    geom_point(
      data = annotated_points_df,
      mapping = aes(
        x = 100 * a1freq,
        y = se
      ),
      col = "grey20",
      size = 0.8
    ) +
    scale_x_continuous(
      name = "Minor Allele Frequency [%]",
      limits = c(0, 50),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "Standard Error Estimate",
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(se_maf_plot)
  
}


# Libraries

library(conflicted)
library(janitor)
library(stringr)
library(glue)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scico)
library(grid)


# Namespace conflicts

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Command line input

args <- commandArgs(TRUE)
 
regenie_output_path <- args[1]
results_folder <- args[2]
docs_folder <- args[3]
analysis_name <- args[4]
p_value_threshold <- args[5]
bp_limit <- args[6]


# DEBUG test data

if (F) {
  
  regenie_output_path <- "/mnt/scratch/marc/moba/2025.12.17/weight_retention/regenie/step_2/bmi_beginning/step2_pop_children_pheno_bmi_diff_end_beginning_chr_chromosome_.regenie.gz"
  results_folder <- "2025.12.17/weight_retention/results/bmi_beginning"
  docs_folder <- "gwas/docs/2025.12.17"
  analysis_name <- "debug_mh"
  p_value_threshold <- 5e-8
  bp_limit <- 500000
  
}

# Housekeeping

doc_folder <- file.path(docs_folder, analysis_name)

md_file <- file.path(doc_folder, glue("pop_{population}_pheno_{pheno}.md"))
figures_folder <- file.path(doc_folder, "figures")
annotation_folder <- file.path(doc_folder, "annotation")

if (!dir.exists(results_folder)) {
  
  dir.create(results_folder, recursive = T)
  
}

if (!dir.exists(doc_folder)) {
  
  dir.create(doc_folder, recursive = T)
  
}

if (!dir.exists(figures_folder)) {
  
  dir.create(figures_folder, recursive = T)
  
}

if (!dir.exists(annotation_folder)) {
  
  dir.create(annotation_folder, recursive = T)
  
}


# GWAS parameters

split <- strsplit(regenie_output_path, "/")[[1]]
split <- strsplit(split[length(split)], "_pop_")[[1]]
split <- strsplit(split[2], "_pheno_")[[1]]

population <- split[1]

split <- strsplit(split[2], "_chr_")[[1]]

pheno <- split[1]

print(paste0(Sys.time(), " - Merging, formatting, and documenting GWAS results for ", analysis_name))

# Plotting parameters

mh_color_1 <- "grey20"
mh_color_2 <- "grey40"

theme_set(theme_bw(base_size = 24))


# Chromosome lengths in GRCh37.p13 (hg19) from Ensembl

chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
chromosome_length <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)
genome_length <- sum(chromosome_length)
chromosome_start <- cumsum(chromosome_length) - chromosome_length
chromosome_middle <- chromosome_start + chromosome_length / 2


# Load the data

association_df <- list()
annotation_df <- list()

print(paste0(Sys.time(), "     Loading regenie results from ", regenie_output_path))

for (chromosome in chromosomes) {
  
  regenie_output_file <- str_replace_all(
    string = regenie_output_path,
    pattern = "_chromosome_", 
    replacement = chromosome
  )
  
  if (file.exists(regenie_output_file)) {
  
  regenie_output <- read.table(
    file = regenie_output_file,
    header = T,
    sep = " ",
    stringsAsFactors = F
  ) %>% 
    clean_names() %>%
    filter(
      !is.na(log10p)
    )
  
  association_df[[length(association_df) + 1]] <- regenie_output
  annotation_df[[length(annotation_df) + 1]] <- get_top_hits(regenie_output)
  
  }
}

association_df <- do.call(rbind, association_df)
annotation_df <- do.call(rbind, annotation_df)

print(paste0(Sys.time(), "     Exporting merged results to ", results_folder))

if (!dir.exists(results_folder)) {

  dir.create(results_folder, recursive = T)

}

write.table(
  x = association_df,
  file = gzfile(file.path(results_folder, glue("pop_{population}_pheno_{pheno}.regenie.gz"))),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  x = annotation_df,
  file = gzfile(file.path(results_folder, glue("pop_{population}_pheno_{pheno}.top_hits.gz"))),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)


# Write documentation

write(
  x = glue("## {pheno} in {population}"),
  file = md_file,
  append = F
)

write(
  x = glue("Association results by regenie for {pheno} in {population}, followed by simple pruning of the hits passing p < {p_value_threshold}.\n"),
  file = md_file,
  append = T
)

if (nrow(regenie_output) == 0) {
  
  write(
    x = glue("**Warning:*** The association results contain no SNP with finite p-value. Please check the log of regenie and the number of cases vs controls."),
    file = md_file,
    append = T
  )
  
  return()
  
}

write(
  x = glue("### Manhattan"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("pop_{population}_pheno_{pheno}_mh.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("pop_{population}_pheno_{pheno}_mh.png"))

write(
  x = glue("![]({relative_figure_path})"),
  file = md_file,
  append = T
)

mh_plot <- get_mh(
  association_df = association_df,
  annotation_df = annotation_df
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(mh_plot)
dummy <- dev.off()

annotation_df_hla <- annotation_df %>% 
  filter(
    chrom == 6 & genpos > 27000000 & genpos < 34000000
  )

annotation_df_not_hla <- annotation_df %>% 
  filter(
    chrom != 6 | genpos <= 27000000 | genpos >= 34000000
  )

write(
  x = glue("### Top hits"),
  file = md_file,
  append = T
)

if (nrow(annotation_df_hla) > 0) {
  
  write(
    x = "> Note: HLA region not included (chr 6, 27-34 Mb)\n", 
    file = md_file, 
    append = T
  )
  
}


write(
  x = paste0("| SNP | chr | bp | allele 0 | allele 1 | allele 1 freq | beta | se | log10p | n |"),
  file = md_file, 
  append = T
)
write(
  x = paste0("| --- | --- | -- | -------- | -------- | ------------- | ---- | -- | ------ | - |"),
  file = md_file, 
  append = T
)

if (nrow(annotation_df_not_hla) > 0) {
  
  for (i in 1:nrow(annotation_df_not_hla)) {
    
    snp <- annotation_df_not_hla$id[i]
    chr <- annotation_df_not_hla$chrom[i]
    bp <- annotation_df_not_hla$genpos[i]
    allele0 <- annotation_df_not_hla$allele0[i]
    allele1 <- annotation_df_not_hla$allele1[i]
    allele1_freq <- annotation_df_not_hla$a1freq[i]
    beta <- annotation_df_not_hla$beta[i]
    se <- annotation_df_not_hla$se[i]
    p <- annotation_df_not_hla$log10p[i]
    n <- annotation_df_not_hla$n[i]
    
    write(
      x = paste0("| ", snp, " | ", chr, " | ", bp, " | ", allele0, " | ", allele1, " | ", allele1_freq, " | ", beta, " | ", se, " | ", p, " | ", n, " |"),
      file = md_file, 
      append = T
    )
    
  }
}

if (nrow(annotation_df_hla) > 0) {
  
  write(
    x = glue("### HLA top hits"),
    file = md_file,
    append = T
  )
  
  write(
    x = paste0("| SNP | chr | bp | allele 0 | allele 1 | allele 1 freq | beta | se | p | n |"),
    file = md_file, 
    append = T
  )
  write(
    x = paste0("| --- | --- | -- | -------- | -------- | ------------- | ---- | -- | - | - |"),
    file = md_file, 
    append = T
  )
  
  
  for (i in 1:nrow(annotation_df_hla)) {
    
    snp <- annotation_df_hla$id[i]
    chr <- annotation_df_hla$chrom[i]
    bp <- annotation_df_hla$genpos[i]
    allele0 <- annotation_df_hla$allele0[i]
    allele1 <- annotation_df_hla$allele1[i]
    allele1_freq <- annotation_df_hla$a1freq[i]
    beta <- annotation_df_hla$beta[i]
    se <- annotation_df_hla$se[i]
    p <- annotation_df_hla$log10p[i]
    n <- annotation_df_hla$n[i]
    
    write(
      x = paste0("| ", snp, " | ", chr, " | ", bp, " | ", allele0, " | ", allele1, " | ", allele1_freq, " | ", beta, " | ", se, " | ", p, " | ", n, " |"),
      file = md_file, 
      append = T
    )
    
  }
}

write(
  x = glue("### Quality Control"),
  file = md_file,
  append = T
)

write(
  x = glue("- QQ plot\n"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("pop_{population}_pheno_{pheno}_qq.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("pop_{population}_pheno_{pheno}_qq.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

qq_plot <- get_qq(
  association_df = association_df
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(qq_plot)
dummy <- dev.off()

write(
  x = glue("- Beta vs. Allele Frequency\n"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("pop_{population}_pheno_{pheno}_beta_af.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("pop_{population}_pheno_{pheno}_beta_af.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

beta_af_plot <- get_beta_maf(
  association_df = association_df,
  annotation_df = annotation_df
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(beta_af_plot)
dummy <- dev.off()

write(
  x = glue("- Standard error vs. Allele Frequency\n"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("pop_{population}_pheno_{pheno}_se_af.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("pop_{population}_pheno_{pheno}_se_af.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

se_af_plot <- get_sd_maf_plot(
  association_df = association_df,
  annotation_df = annotation_df
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(se_af_plot)
dummy <- dev.off()

