# GWAS

## Inputs

### Data

For this workflow, both bed and bgen genotypes should be available. They will likely be the respective outputs of the [Combining Datasets](@ref) and [Genotypes Imputation](@ref) workflows, optionally merged with [Combining with UK Biobank](@ref) depending on your project of interest.

### GWAS Variables

The variables to be used for the GWAS are specified in a YAML file like so.

```yaml
group: ["ANCESTRY"]
phenotype: "SEVERE_COVID_19"
covariates: ["AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE"]
```

where:

- A GWAS will be run independently for each group defined by `group`
- The phenotype is the trait of interest, currently only `SEVERE_COVID_19`
- The covariates are additional covariates that will enter the GWAS to improve precision

## Workflow Parameters

This is the list of all the pipeline's parameters, they can be set in the `run.config` file under the `params` section.

- `GENOTYPES_PREFIX`: Prefix to genotypes in PLINK BED format (likely the output of the [Combining Datasets](@ref) workflow).
- `IMPUTED_GENOTYPES_PREFIX` = Prefix to imputed genotypes in PGEN format (likely the output of the [Genotypes Imputation](@ref) workflow).
- `COVARIATES`: Path to covariate file (likely the output of the [Combining Datasets](@ref) workflow)
- `N_PCS (default: 10)`: Number of principal components to compute.
- `PCA_APPROX` (default: true): Whether PCA is performed via approximation [see](https://www.cog-genomics.org/plink/2.0/strat)
- `MIN_GROUP_SIZE (default: 100)`: Minimum number of samples in a group to proceed to effect size estimation.
- `VARIABLES_CONFIG (default: assets/variables.yaml)`: File containing the declaration of groups, phenotypes and covariates for the GWAS (see [GWAS Variables](@ref)).
- `REGENIE_MAF (default: 0.01)`: Minor allele frequency for a variant to enter the GWAS.
- `REGENIE_MAC (default: 10)`: Minor allele count for a variant to enter the GWAS.
- `REGENIE_BSIZE (default: 1000)`: Regenie's block size (see the [regenie docs](https://rgcgithub.github.io/regenie/))

## Running The Workflow

If the previous steps have been completed successfully you can run:

```bash
./run.sh GWAS
```

## Outputs

All outputs are produced in `PUBLISH_DIR`, the main outputs of the workflow are:


## Workflow DAG

```@raw html
<iframe src="../assets/gwas_dag.html" width="100%" height="800px"></iframe>
```