# Combining GenOMICC Datasets

This workflow combines the various genetic data sources available into a unified dataset.

## Inputs

Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This section describes the input data, the corresponding workflow parameters are described below.

### Genetic Data

The pipeline requires genotyping arrays and optionally whole-genome sequencing (wgs) data.

#### Genotyping Arrays

There are three filesets and three corresponding subfolders:

- The r8 release: Corresponds to genotyping data generated before 2021. The genotyping chip was the Illumina GSA-MD-24v3-0_A1 and the associated genome build GRCh37. The corresponding subfolder is usually named `wp5-gwas-r8-under90excl_2021Sep16`.
- The 2021-2023 release:  Corresponds to genotyping data generated between 2021 and 2023. It was also genotyped using the Illumina GSA-MD-24v3-0_A1 chip and the genome build is also GRCh37. The corresponding subfolder is usually named `20210920_20231206_QC_VFinal`.
- The 2024-now release: Corresponds to the latest fileset. The Illumina GSA-MD-48v4-0_A1 chip was used and the genome build is GRCh38. The corresponding subfolder is usually named `20240604_20240610_QC_VFinal`.

In each subfolder, there is a plink subfolder which contains the actual genotypes, the following tables summarises the above

| Brief Description | Period Begin | Period End | Genotyping Array | Genome Build | Directory | Genotypes Prefix |
| :--------: | :------------: | :----------: | :----------------: | :------------: | :------------: | :------------: |
| Prehistoric r8 release | 04/05/2020   | 30/08/2021 | GSA-MD-24v3-0_A1               | GRC37        | wp5-gwas-r8-under90excl_2021Sep16 | PLINK_190921_0906/wp5-gwas-r8-under90excl_2021Sep16 |
| Before 2024 microarray | 20/09/2021   | 06/12/2023 | GSA-MD-24v3-0_A1 | GRC37 | 20210920_20231206_QC_VFinal | PLINK_040724_0954/20210920_20231206_QC_VFinal |
| Since 2024 microarray  | 04/06/2024   | 10/06/2024 | GSA-MD-48v4-0_A1 | GRC38 | 20240604_20240610_QC_VFinal | PLINK_040724_0114/20240604_20240610_QC_VFinal |

So, in the example above, the `R8_GENOTYPES` workflow parameter (see below) should point to `wp5-gwas-r8-under90excl_2021Sep16/PLINK_190921_0906/wp5-gwas-r8-under90excl_2021Sep16`

#### Whole Genome Sequencing

The wgs GVCF files are all located in a `wgs-reheadered` folder.

### External Resources

As well as our in-house data, the pipeline depends on external reference data. In principle these files should already be present on ODAP and there is nothing you need to do.

#### The 1000 GP

- All VCF files and indexes present in [this FTP folder](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
- The associated 1000 GP [pedigree file](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt).

These files should be stored in the same folder which is defined by the `KGP_DIR (default: /mnt/odap-beegfs/software/gwas-resources/1000-genomes-HC)` Nextflow parameter.

#### GATK

- The [reference genome](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta) published by the Broad Institute.

This should be in a folder defined by the `GATK_DIR (default: /mnt/odap-beegfs/software/gwas-resources/gatk)` Nextflow parameter.
  
## Running The Workflow

If the previous steps have been completed successfully you can run:

```bash
taskset -c 998 nextflow run main.nf -entry CombineGeneticDatasets -c run.config -profile odap -resume -with-report -with-trace
```

## Outputs

All outputs are produced in `PUBLISH_DIR` (defaults to `results`), the main outputs of the workflow are:

- `report.md`: A report of the pipeline execution
- `genotypes.aggregated.qced.final.{bed,bim,fam}`: The aggregated genotypes.
- `covariates.inferred.csv`: The covariates inferred from the genotypes (ancestry, PCs).
- `ancestry/kgp_shared/kgp.merged.unrelated.shared.{bed,bim,fam}`: The 1000 Genome Project genotypes corresponding to the variants in the GenOMICC cohort (can be used later on in [Merging the GenOMICC and UK Biobank Cohorts](https://github.com/baillielab/ukb-genomicc-workflows))

## Workflow Parameters

This is the list of all the pipeline's parameters, they can be set in the `run.config` file under the `params` section.

### Input Files

These are project specific and need to be provided:

- `R8_GENOTYPES`: Prefix to release r8 genotypes (before 2021).
- `BEFORE_2024_GENOTYPES`: Prefix to genotypes released between 2021-2023.
- `SINCE_2024_GENOTYPES`: Prefix to genotypes released after 2024.
- `WGS_GVCFS` (optional): Prefix to whole genome sequencing files.

### External Inputs Parameters

These are already set if you are using the `odap` profile.

- `RESOURCES_DIR` (default: ./assets/resources"): Path to all external resources.
- `KGP_DIR`: Path to the 1000 Genome Project specific resources (see [The 1000 GP](@ref)).
- `GATK_DIR`: Path to GATK specific resources (see [GATK](@ref)).
- `GRC37_TO_GRC38_CHAIN_FILE` (default: "./assets/hg19ToHg38.over.chain.gz"): Path to chain file used to liftover the GRCh37 genotypes to GRCh38.

### QC Parameters

- `QC_GENOTYPE_MISSING_RATE` (default: 0.02): Maximum missing rate per variant across all individuals. Variants above the threshold are dropped.
- `QC_INDIVIDUAL_MISSING_RATE` (default: 0.02): Maximum missing rate per individual across genotypes. Individuals above the threshold are dropped.
- `QC_HWE_P` (default: 1e-10): Used to identify potential technical artifacts and drop variants.
- `QC_HWE_K` (default: 0.001): Used together with `QC_HWE_P`
- `PCA_APPROX` (default: true): Whether PCA is performed via approximation [see](https://www.cog-genomics.org/plink/2.0/strat)
- `FILTER_HIGH_LOADINGS_VARIANTS` (default: false): Whether to drop variants with high PCA loadings. If the loadings plots exhibits a high peak you may want to turn that on.
- `ANCESTRY_THRESHOLD` (default: 0.5): For each individual, the most likely ancestry estimate should be greater than this threshold otherwise the individual is marked as admixed.
- `ANCESTRY_PROGRAM` (default: scope): Program to use for ancestry estimation. Either `scope` or `admixture`.

### Output Directories Parameters

- `PUBLISH_DIR` (default: "results"): Top level directory where outputs will be output.
- `KGP_PUBLISH_DIR` (default: "results/kgp"): Where 1000 Genome Project data will be output.
- `ARRAY_GENOTYPES_PUBLISH_DIR` (default: "results/array_genotypes"): Where data associated with the genotyping arrays will be output.
- `WGS_PUBLISH_DIR` (default: "results/wgs"): Where data associated with the whole-genome sequencing data will be output.
- `GATK_PUBLISH_DIR` (default: "results/gatk"): Where data associated with GATK requirements will be output.
- `MERGED_PUBLISH_DIR` (default: "results/merged"): Where the merged genetic data will be output.

## Current Limitations

These are current limitations of the aggregation workflow:

- Only chromosomes 1 to 22 are processed.
- Only bi-allelic SNPs are used.

## Workflow DAG

```@raw html
<iframe src="../assets/combining_datasets_dag.html" width="100%" height="2000px"></iframe>
```