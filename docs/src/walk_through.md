# Walk Through

In this section we discuss how to regenerate results from an analysis from scratch. Each step consists in a separate workflow for which a dedicated documentation is referenced. While workflows are platform agnostic, due to complex data permissions and governance, each workflow is typically run on a dedicated platform.

## Step 1: Aggregating Genetic Datasets

- workflow: [Combining GenOMICC Datasets](@ref)
- platform: [ODAP](@ref)

## Step 2: Imputation of Genotypes

- workflow: [TOPMed Imputation](https://github.com/olivierlabayle/nf-topmed-imputation)
- platform: [Eddie](https://information-services.ed.ac.uk/research-support/research-computing/ecdf/high-performance-computing)

## Step 3: Merging GenOMICC and the UK Biobank

- workflow: [Merging the GenOMICC and UK Biobank Cohorts](https://github.com/baillielab/ukb-genomicc-workflows)
- platform: [UKB RAP](https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/)

## Step 4: GWAS

- workflow: [WDL-GWAS](https://github.com/olivierlabayle/WDL-GWAS)
- platform: [UKB RAP](https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/)
