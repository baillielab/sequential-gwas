# Combining Datasets

This workflow combines the various available data sources into a unified dataset.

## Inputs

Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This page aims at making it clear what are the inputs and outputs to the pipeline.

All file paths are currently reported with respect to my a015 project but this will be standardized in the future.

### External Resources

As well as our in-house data, the pipeline depends on external reference data. In principle these files should already be present on ODAP and there is nothing you need to do.

#### The 1000 GP

- All VCF files and indexes present in [this FTP folder](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
- The associated 1000 GP [pedigree file](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt).

These files should be stored in the same folder which is defined by the `KGP_DIR (default: /mnt/odap-beegfs/software/gwas-resources/1000-genomes-HC)` Nextflow parameter.

#### GATK

- The [reference genome](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta) published by the Broad Institute.

This should be in a folder defined by the `GATK_DIR (default: /mnt/odap-beegfs/software/gwas-resources/gatk)` Nextflow parameter.

### Genetic Data

The pipeline requires both genotyping arrays and whole-genome sequencing (wgs) data. These should be all be located within a `pre-qc` directory and organised as follows:

#### Genotyping Arrays

There are three filesets and three corresponding subfolders:

- The r8 release: Corresponds to genotyping data generated before 2021. The genotyping chip was the Illumina GSA-MD-24v3-0_A1 and the associated genome build GRCh37. The corresponding subfolder is `wp5-gwas-r8-under90excl_2021Sep16`.
- The 2021-2023 release:  Corresponds to genotyping data generated between 2021 and 2023. It was also genotyped using the Illumina GSA-MD-24v3-0_A1 chip and the genome build is also GRCh37. The corresponding subfolder is `20210920_20231206_QC_VFinal`.
- The 2024-now release: Corresponds to the latest fileset. The Illumina GSA-MD-48v4-0_A1 chip was used and the genome build is GRCh38. The corresponding subfolder is `20240604_20240610_QC_VFinal`.

The following tables summarises the above

| Brief Description | Period Begin | Period End | Genotyping Array | Genome Build | Directory | Genotypes File |
| :--------: | :------------: | :----------: | :----------------: | :------------: | :------------: | :------------: |
| Prehistoric r8 release | 04/05/2020   | 30/08/2021 | GSA-MD-24v3-0_A1               | GRC37        | wp5-gwas-r8-under90excl_2021Sep16 | PLINK_190921_0906/wp5-gwas-r8-under90excl_2021Sep16.ped |
| Before 2024 microarray | 20/09/2021   | 06/12/2023 | GSA-MD-24v3-0_A1 | GRC37 | 20210920_20231206_QC_VFinal | PLINK_040724_0954/20210920_20231206_QC_VFinal.ped |
| Since 2024 microarray  | 04/06/2024   | 10/06/2024 | GSA-MD-48v4-0_A1 | GRC38 | 20240604_20240610_QC_VFinal | PLINK_040724_0114/20240604_20240610_QC_VFinal.ped |

#### Whole Genome Sequencing

The wgs GVCF files are all located in a `wgs-reheadered` folder.

### Covariates Data

The covariate file is a CSV file in the `pre-qc` directory and named `covariates.csv` with the following columns:

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
  
## Running The Workflow

If the previous steps have been completed successfully you can run:

```bash
bash run.sh
```

## Pipeline parameters

This is the list of all the pipeline's parameters. In principle they don't need to be changed if the conventions in this documentation have been respected and are up to date. Otherwise, please feel free to open an issue.

### Output Directories Parameters

- `PUBLISH_DIR` (default: "results"): Top level directory where outputs will be output.
- `KGP_PUBLISH_DIR` (default: "results/kgp"): Where 1000 Genome Project data will be output.
- `ARRAY_GENOTYPES_PUBLISH_DIR` (default: "results/array_genotypes"): Where data associated with the genotyping arrays will be output.
- `WGS_PUBLISH_DIR` (default: "results/wgs"): Where data associated with the whole-genome sequencing data will be output.
- `GATK_PUBLISH_DIR` (default: "results/gatk"): Where data associated with GATK requirements will be output.
- `MERGED_PUBLISH_DIR` (default: "results/merged"): Where the merged genetic data will be output.

### External Inputs Parameters

- `RESOURCES_DIR` (default: ./assets/resources"): Path to all external resources.
- `KGP_DIR` (default: "./assets/kgp"): Path to the 1000 Genome Project specific resources.
- `GATK_DIR` (default: "./assets/gatk"): Path to GATK specific resources.
- `GRC37_TO_GRC38_CHAIN_FILE` (default: "./assets/hg19ToHg38.over.chain.gz"): Path to chain file used to liftover the GRCh37 genotypes to GRCh38.

### Basic QC Parameters

- `QC_GENOTYPE_MISSING_RATE` (default: 0.1): Maximum missing rate per variant across all individuals. Variants above the threshold are dropped.
- `QC_INDIVIDUAL_MISSING_RATE` (default: 0.1): Maximum missing rate per individual across genotypes. Individuals above the threshold are dropped.
- `QC_HWE_P` (default: 1e-5): Used to identify potential technical artifacts and drop variants.
- `QC_HWE_K` (default: 0.001): Used together with QC_HWE_P

## DAG

```@raw html
<iframe src="./assets/combining_datasets_dag.html" width="100%" height="800px"></iframe>
```