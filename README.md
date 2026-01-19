# GWAS base pipeline

This is an example pipeline to run GWASes on the MoBa genotypes. It is based on the HDGB release of the genotypes organized on the `bergomo` lab on HUNT cloud.

## Requirements
The pipeline should be used with phenotypes with version posterior to `25-12-12`.

## Content
1. `phenotypes`: This folder contains a quatro document for the processing and documentation of phenotypes for your GWAS. It creates a table that is picked by the GWAS pipeline.
2. `gwas`: This folder contains a snakemake file that runs the GWAS on the table of phenoptypes produced in the previous step.

## Instructions
- Step 0: Fork this repository or copy its content to your own repository. Pull the repository on HUNT cloud via ssh.
- Step 1: Open the quatro document `phenotypes/phenotypes.qmd` using RStudio in workbench. Extract the summary statistics needed for your study as well as plots documenting the distribution of values. Note the path to the table of phenotypes.
- Step 2: Edit the GWAS pipeline settings in `gwas/analysis.yaml` and execute the file `gwas/gwas_pipeline_hdgb.snakemake` via ssh.

## Known issues
- The pruning of hits currently does not take LD into account, this will be implemented at a later stage.
- Regenie is used only one phenotype at a time to allow using different sets of covariates for each phenotype. Grouping phenotypes would allow speed improvements but requires structural changes to the pipeline.
- GWAS catalog results are currently not available.

## Example
The repository contains results from an example on birth weight.
- [phenotypes.qmd](phenotypes/phenotypes.qmd): contains the phenotypes handling of the birth weight analysis.
- [children](gwas/docs/2026.01.06/birth_weight/pop_children_pheno_weight_birth.md): GWAS against the genotypes of the children.
- [mothers](gwas/docs/2026.01.06/birth_weight/pop_mothers_pheno_weight_birth.md): GWAS against the genotypes of the mothers.
- [fathers](gwas/docs/2026.01.06/birth_weight/pop_fathers_pheno_weight_birth.md): GWAS against the genotypes of the fathers.
- [parents](gwas/docs/2026.01.06/birth_weight/pop_parents_pheno_weight_birth.md): GWAS against the genotypes of the parents.


