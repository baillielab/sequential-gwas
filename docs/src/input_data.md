# Input Data

Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This page aims at making it clear what are the inputs and outputs to the pipeline.

All file paths are currently reported with respect to my a015 project but this will be standardized in the future.

## Genetic Data

Depending on the time period, individuals may have been genotyped, whole genome sequenced or both.

### Genotyping Arrays

| Brief Description | Period Begin | Period End | Genotyping Array | Genome Build | Directory | Genotypes File |
| :--------: | :------------: | :----------: | :----------------: | :------------: | :------------: | :------------: |
| Prehistoric r8 release | 04/05/2020   | 30/08/2021 | GSA-MD-24v3-0_A1               | GRC37        | wp5-gwas-r8-under90excl_2021Sep16 | PLINK_190921_0906/wp5-gwas-r8-under90excl_2021Sep16.ped |
| Before 2024 microarray | 20/09/2021   | 06/12/2023 | GSA-MD-24v3-0_A1 | GRC37 | 20210920_20231206_QC_VFinal | PLINK_040724_0954/20210920_20231206_QC_VFinal.ped |
| Since 2024 microarray  | 04/06/2024   | 10/06/2024 | GSA-MD-48v4-0_A1 | GRC38 | 20240604_20240610_QC_VFinal | PLINK_040724_0114/20240604_20240610_QC_VFinal.ped |

The Illumina manifest files corresponding to each array were downloaded from [Illumina's website](https://emea.support.illumina.com/array/array_kits/infinium-global-screening-array/downloads.html) and are stored in `assets`. The description of the manifest columns can be found [here](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565).

### Whole Genome Sequencing

Currently the data available to me is a subset of ~500 000 SNPs split in files of groups of individuals in `/odp-beegfs/a015/linked_data/preqc/gws-genotyped`.

### Erola's Historical Data (This paragraph will be deleted)

As currently understood and relative to the data filesystem in `/odp-beegfs/a015/postqc/array/`:

- `imputation-output/chr*.info.gz` are outputs from TopMed imputation services from `erola_scripts/downloadTOPMED.sh`. These originate from merged/liftover genotyped files that have been lost.
- `imputation-output/Genomicc_chr*.txt` are probable outputs from `erola_scripts/filtervcf.sh` based on the above even though there seems to be an intermediate vcf file that was deleted.
- In `bgen/` are the likely output of `erola_scripts/vcftobgen.sh`.

## Covariates Data

The covariate file is a CSV file in `/odp-beegfs/a015/linked_data/postqc/array/a015_covariates.csv` with the following columns:

- ODAP_ID: Unique invidual identifier (same as the IID/FID in the genetic data).
- COHORTS: List of all cohorts the individual belongs to.
- AGE_AT_RECRUITMENT: Age of the individual at recruitment. If an individual was recruited in two cohorts, typically (GenOMICC and ISARIC), GenOMICC takes precedence.
- REPORTED_SEX: Sex reported by an individual.
- PRIMARY_DIAGNOSIS: Only for GenOMICC patients and obtained via REDCap.
- DIAGNOSES: List of diagnoses inferred from REDCap
- ISARIC_SEVERITY_SCORE: Only for ISARIC individuals, severity of infection.
- GENETIC_SAMPLE_TYPE: Blood or Saliva.
- CALL_RATE: Genotyping call rate (saliva samples are believed to have lower quality).
- GENETIC_MEASUREMENT_TECHNOLOGY: At this time either of (GSA-MD-24v3-0_A1, GSA-MD-48v4-0_A1, WGS).
- DATE_AT_RECRUITMENT: Time marker of individual recruitement.