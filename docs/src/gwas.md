# GWAS

This WDL workflow performs a genome-wide association study across all defined subgroups and phenotypes. It can either be run on [ODAP](@ref) or on the [UKB RAP](@ref). 

## Inputs

Depending on your phenotypes of interest, some of the input files to the workflow will be different. For a mortality study, typically only GenOMICC is used and the workflow can be run on ODAP. For a susceptibility study, typically both GenOMICC and UK Biobank participants are used and the workflow will be run on the UK Biobank RAP. As such some inputs to the workflow will either be outputs of [Combining GenOMICC Datasets](@ref) or [Merging the GenOMICC and UK Biobank Cohorts](@ref).

- `docker_image` (default: olivierlabayle/genomicc:scope): The docker image used to run the workflow.
- `covariates_file`: TODO
- `genotypes`: 
  - For a mortality study, these are the `genotypes.aggregated.qced.final.{bed,bim,fam}` output by [Combining GenOMICC Datasets](@ref).
  - For a susceptibility study, these are the genotypes from [Merging Cohorts Outputs](@ref).
- `imputed_genotypes`:
  - For a mortality study, these are the imputed genotypes from [Imputation Outputs](@ref).
  - For a susceptibility study, these are the imputed genotypes from [Merging Cohorts Outputs](@ref).
- `groupby`: A set of variables used to stratify individuals for which a GWAS will be run independently. If empty, the full dataset is used.
- `covariates`: A set of covariates used to adjust for confounding or increase power in the association testing step. Product of variables can be defined using the `_x_` syntax, for example: ["AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE"].
- `phenotypes`: The set of binary phenotypes for which a GWAS will be run independently.
- `min_group_size` (default: 100): When using `groupby`, groups must consist of at least `min_group_size` individuals.
- `high_ld_regions`: File containing high LD regions to be excluded when performing LD pruning for PCA. The file is stored in `assets/exclude_b38.txt` and needs to be uploaded to the RAP.
- `ip_values` (default: "1000 50 0.05"): Values used to create independent genotypes for PCA (see [here](https://www.cog-genomics.org/plink/2.0/ld)).
- `npcs` (default 10): Number of principal components to use to account for population structure.
- `approx_pca` (default: true): Whether to use an approximation to the PCA algorithm (see [here](https://www.cog-genomics.org/plink/2.0/strat)).
- `maf` (default: 0.01): Minor allele frequency used to filter variants entering the GWAS.
- `mac` (default: 10): Minor allele count used to filter variants entering the GWAS.
- `regenie_cv_folds` (default: 5): Number of folds for Regenie step 1.
- `regenie_bsize` (default 1000): Regenie block size.

For Regenie's options see the [online documentation](https://rgcgithub.github.io/regenie/options/).

The inputs need to be filled within the `rap_workflows/gwas/inputs.json` file.

## Running On UKB RAP

First you need to compile the WDL workflow and upload it to the RAP, this can be done with the following:

```bash
export DX_COMPILER_PATH=/Users/olabayle/dxCompiler/dxCompiler-2.13.0.jar
export PROJECT_ID=project-J0pkqyQJpYQ133JG1p2J1qzv
java -jar $DX_COMPILER_PATH compile rap_workflows/gwas/workflow.wdl -f -project $PROJECT_ID -folder /workflows/gwas -inputs rap_workflows/gwas/inputs.json
```

where the `DX_COMPILER_PATH` and `PROJECT_ID` have to be set appropriately. The compiler might output some warnings like `missing input for non-optional parameter` but you can ignore these.

!!! warning "Updating workflows"
    At this point in time it seems like compiling multiple times the same workflow does not replace the old files. You will need to manually erase them from the RAP.

Then, you can run the workflow with the following command

```bash
dx run -y \
-f rap_workflows/gwas/inputs.dx.json \
--priority high \
--destination /gwas_outputs/ \
/workflows/gwas/gwas
```

## Outputs

TODO