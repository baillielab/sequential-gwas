# Input Data

Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This page aims at making it clear what are the inputs and outputs to the pipeline.

All file paths are currently reported with respect to my a015 project but this will be standardized in the future.

## External Resources

As well as our in-house data, the pipeline depends on external reference data. In principle these files should already be present on ODAP and there is nothing you need to do.

### The 1000 GP

- All VCF files and indexes present in [this FTP folder](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
- The associated 1000 GP [pedigree file](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt).

These files should be stored in the same folder which is defined by the `KGP_DIR (default: /mnt/odap-beegfs/software/gwas-resources/1000-genomes-HC)` Nextflow parameter.

### GATK

- The [reference genome](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta) published by the Broad Institute.

This should be in a folder defined by the `GATK_DIR (default: /mnt/odap-beegfs/software/gwas-resources/gatk)` Nextflow parameter.

## Genetic Data

The pipeline requires both genotyping arrays and whole-genome sequencing (wgs) data. These should be all be located within a `pre-qc` directory and organised as follows:

### Genotyping Arrays

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



### Whole Genome Sequencing

The wgs GVCF files are all located in a `wgs-reheadered` folder.

## Covariates Data

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