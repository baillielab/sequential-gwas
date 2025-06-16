# Merging the GenOMICC and UK Biobank Cohorts

This workflow yields combined:

- genotypes in PLINK bed format: obtained by subsampling the TOPMed UKB imputed genotypes to keep the GenOMICC variants.
- imputed genotypes in PLINK pgen format: ontained by merging UKB and GenOMICC TOPMed UKB imputed genotypes.
- covariates containing the current age and sex of individuals
- ancestry estimates: via the 1000 Genome Project

## Inputs

This section describes the main inputs to the workflow, see further down for the associated Nextflow parameters.

### GenOMICC Data

For this workflow, genotypes, imputed genotypes are required. They will be the respective outputs of the [Combining Datasets](@ref) and [Genotypes Imputation](@ref) workflows.

You will also need two covariates files:

- The `COVARIATES` file contains the usual covariates information including phenotypes and is typically provided by Dominique on a per project basis. Unfortunately, at this point, the specifications for this file have not been clearly defined. The columns of the file I have received are: `genotype_file_id`, `age_years`, `sex`, `case_or_control`, `cohort`, `severe_cohort_primary_diagnosis`, `isaric_cohort_max_severity_score`.
- An optional `INFERRED_COVARIATES` contains covariates information that is inferred from the genetic datasets during the [Combining Datasets](@ref) workflow. This is the case for ancestry estimates for instance.

### 1000 Genome Project

You will also need the 1000 Ggenome Project genotypes in plink format, this is typically output by the [Combining Datasets](@ref) workflow, but can also be obtained by running the yet undocumented Nextflow `KGP` workflow in this repository.

### UK Biobank data

The UK-Biobank data is on the RAP, there is nothing you need to do but to setup your account and follow the instructions in the [Working with DNANexus](@ref) section.

## Workflow Parameters

In `wdl/ukb_merge/inputs.json`, find and replace all instances of the project-ID with your own.

## Running The Workflow

First you need to compile the WDL workflow and upload it to the RAP, this can be done by the following:

```bash
export DX_COMPILER_PATH=/Users/olabayle/dxCompiler/dxCompiler-2.13.0.jar
export PROJECT_ID=project-J0pkqyQJpYQ133JG1p2J1qzv
java -jar $DX_COMPILER_PATH compile wdl/ukb_merge/workflow.wdl -f -project $PROJECT_ID -folder /workflows -inputs wdl/ukb_merge/inputs.json
```

where the `DX_COMPILER_PATH` and `PROJECT_ID` have to be set appropriately. The compiler might output some warnings like `missing input for non-optional parameter` but you can ignore these.

Then, you can run the workflow with the following command

```bash
dx run /workflows/merge_ukb_and_genomicc -f wdl/ukb_merge/inputs.dx.json
```

## Outputs
