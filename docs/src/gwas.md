# GWAS


## Inputs

This section describes the main inputs to the workflow, see further down for the associated Nextflow parameters.

### Genetic Data

For this workflow, both genotypes and imputed genotypes are required. They will likely be the respective outputs of the [Combining Datasets](@ref) and [Genotypes Imputation](@ref) workflows, optionally merged with [Combining with UK Biobank](@ref) depending on your project of interest.

### Covariates

There are two main types of covariates files you can use here. 

- The `COVARIATES` file contains the usual covariates information including phenotypes and is typically provided by Dominique on a per project basis. Unfortunately, at this point, the specifications for this file have not been clearly defined. The columns of the file I have received are: `genotype_file_id`, `age_years`, `sex`, `case_or_control`, `cohort`, `severe_cohort_primary_diagnosis`, `isaric_cohort_max_severity_score`.
- An optional `INFERRED_COVARIATES` contains covariates information that is inferred from the genetic datasets during the [Combining Datasets](@ref) workflow. This is the case for ancestry estimates for instance.

### GWAS Variables

You can define the phenotypes, covariates and groups via an external YAML file which looks like the following. 

```yaml
groupby: ["ANCESTRY_ESTIMATE"]
phenotypes: ["SEVERE_COVID_19"]
covariates: ["AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE"]
```

In this example:

- An independent GWAS will be run for each ancestry group in the `ANCESTRY_ESTIMATE` column. 
- There is only one phenotype: `SEVERE_COVID_19`. 
- The covariates to are ``AGE``, ``SEX``, ``AGE^2`` and ``AGE \cdot AGE``.

#### Groups

The `groupby` variables defines groups for which an independent GWAS will be run. If the ``groupby`` key is not present in the YAML file, only one GWAS will be run for the whole population.

#### Phenotypes

The `phenotype` variables must either be `case_or_control` or have an existing definition based on other columns in the `COVARIATES` file (hardcoded), at the moment only `SEVERE_COVID_19` has been defined.

#### Covariates 

The `covariates` variables are used as adjustment variables in the regression model. In REGENIE step 2, principal components are also added in a leave-one-chromosome-out scheme to the covariates defined here. 

Again, due to the ill-definition of the covariates file, at the moment, only a subset of non-inferred covariates can be processed: `AGE` and `SEX`. Categorical covariates are one-hot encoded.

!!! note "Interaction Terms"
    Interaction terms can be added to the covariate list with the ``_x_`` delimiter.

## Workflow Parameters

This is the list of all the pipeline's parameters, they can be set in the `run.config` file under the `params` section.

- `GENOTYPES_PREFIX`: Prefix to genotypes in PLINK BED format (likely the output of the [Combining Datasets](@ref) workflow).
- `IMPUTED_GENOTYPES_PREFIX`: Prefix to imputed genotypes in PGEN format (likely the output of the [Genotypes Imputation](@ref) workflow).
- `COVARIATES`: Path to covariate file (likely the output of the [Combining Datasets](@ref) workflow)
- `INFERRED_COVARIATES`: Optional, path to covariates inferred from genetic data during the [Combining Datasets](@ref) workflow.
- `N_PCS (default: 10)`: Number of principal components to compute.
- `PCA_APPROX` (default: true): Whether PCA is performed via approximation [see](https://www.cog-genomics.org/plink/2.0/strat)
- `MIN_GROUP_SIZE (default: 100)`: Minimum number of samples in a group to proceed to effect size estimation.
- `VARIABLES_CONFIG (default: assets/variables.yaml)`: File containing the declaration of groups, phenotypes and covariates for the GWAS (see [GWAS Variables](@ref)).
- `REGENIE_MAF (default: 0.01)`: Minor allele frequency for a variant to enter the GWAS.
- `REGENIE_MAC (default: 10)`: Minor allele count for a variant to enter the GWAS.
- `REGENIE_BSIZE (default: 1000)`: Regenie's block size (see the [regenie docs](https://rgcgithub.github.io/regenie/))

## Running The Workflow

If the previous steps have been completed successfully you can run on your platform (e.g. `PLATFORM=eddie`):

```bash
nextflow run main.nf -entry GWAS -resume -profile PLATFORM -with-trace -with-report -c run.config
```

## Outputs

All outputs are produced in `PUBLISH_DIR (default: results)`, for each group a separate folder contains:

- plots: Manhattan and QQ plots
- results: A CSV file of summary statistics

## Workflow DAG

```@raw html
<iframe src="../assets/gwas_dag.html" width="100%" height="800px"></iframe>
```