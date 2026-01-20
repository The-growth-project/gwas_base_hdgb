
##
#
# Writes the main documentation.
#
##


# Libraries

library(conflicted)
library(yaml)
library(glue)
library(dplyr)
library(ggplot2)
library(ggside)
library(scico)
library(grid)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)

# General parameters
theme_set(theme_bw(base_size = 14))

# Import gwas_settings file

gwas_settings <- read_yaml("gwas/analysis.yaml")


# Write readme

output_folder <- glue("gwas/docs/{gwas_settings$release}")

readme_file <- file.path(output_folder, "readme.md")

write(
  x = glue("# {gwas_settings$name}\n"),
  file = readme_file,
  append = F
)
write(
  x = glue("{gwas_settings$description}\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("The documentation corresponds to the analyses version `{gwas_settings$release}`.\n"),
  file = readme_file,
  append = T
)

for (analysis_id in names(gwas_settings$analyses)) {
  
  # Link to analysis-specific md
  
  analysis <- gwas_settings$analyses[[analysis_id]]
  
  write(
    x = glue("### {analysis$name}\n"),
    file = readme_file,
    append = T
  )
  
  write(
    x = glue("{analysis$description}\n\n"),
    file = readme_file,
    append = T
  )
  
  for (population in analysis$populations) {
    
    relative_path <- glue("{gwas_settings$release}/pop_{population}_pheno_{analysis$phenotype}.md")
    
    write(
      x = glue("- [{population}]({relative_path}): GWAS of {analysis$phenotype_name} against the genome of {population}.\n\n"),
      file = readme_file,
      append = T
    )
  }
}


write(
  x = "#### License\n",
  file = readme_file,
  append = T
)
write(
  x = "Unless otherwise specified in specific files, this work, the associated source code, and results are released under a [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).\n",
  file = readme_file,
  append = T
)
write(
  x = "![CC BY 4.0 License Logo](https://i.creativecommons.org/l/by/4.0/88x31.png)\n\n",
  file = readme_file,
  append = T
)
