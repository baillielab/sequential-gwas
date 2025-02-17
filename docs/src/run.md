# Running The Workflow

## Dependencies

### Platform: ODAP

Since the data is sensitive, it is meant to run on the [ODAP](https://odap.ac.uk/) platform (at least for now). To get access, you will need to contact the relevant person, at the moment [dominique.mccormick@ed.ac.uk](mailto:dominique.mccormick@ed.ac.uk). More information on accessing ODAP can be found [here](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide).

### Data

Once access to ODAP is given you will also need access to the [Input Data](@ref). This is also handled by Dominique at the moment.

### Software

In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:

- Nextflow 24.10.3 ([see how to setup](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/nexflow))
- Singularity 3.9.4 (Should be ready to use)

## Running The Workflow

If the previous steps have been completed successfully you can run:

```bash
nextflow run main.nf -profile odap -resume
```

## Pipeline parameters

We now list all the pipeline parameters. In principle they don't need to be changed if the conventions in this documentation have been respected and are up to date. Otherwise, please feel free to run an issue.

### Output Directories Parameters

- PUBLISH_DIR (default: "results"): Top level directory where outputs will be output.
- KGP_PUBLISH_DIR (default: "${params.PUBLISH_DIR}/kgp"): Where 1000 Genome Project data will be output.
- ARRAY_GENOTYPES_PUBLISH_DIR (default: "${params.PUBLISH_DIR}/array_genotypes"): Where data associated with the genotyping arrays will be output.
- WGS_PUBLISH_DIR (default: "${params.PUBLISH_DIR}/wgs"): Where data associated with the whole-genome sequencing data will be output.
- GATK_PUBLISH_DIR (default: "${params.PUBLISH_DIR}/gatk"): Where data associated with GATK requirements will be output.
- MERGED_PUBLISH_DIR (default: "${params.PUBLISH_DIR}/merged"): Where the merged genetic data will be output.

### External Inputs Parameters

- RESOURCES_DIR (default: ${projectDir}/assets/resources"): Path to all external resources.
- KGP_DIR (default: "${params.RESOURCES_DIR}/kgp"): Path to the 1000 Genome Project specific resources.
- GATK_DIR (default: "${params.RESOURCES_DIR}/gatk"): Path to GATK specific resources.
- GRC37_TO_GRC38_CHAIN_FILE (default: "${projectDir}/assets/hg19ToHg38.over.chain.gz"): Path to chain file used to liftover the GRCh37 genotypes to GRCh38.
- VARIANTS_TO_FLIP_GRC38 (default: "${projectDir}/assets/GSA-48v4-0_20085471_D2-minus-strand.txt"): THis is not used at the moment and probably needs to be dropped.
- VARIANTS_TO_FLIP_GRC37 (default: "${projectDir}/assets/GSA-24v3-0_A1-minus-strand.txt"): same as above.

### Basic QC Parameters

- QC_GENOTYPE_MISSING_RATE (default: 0.1): Maximum missing rate per variant across all individuals. Variants above the threshold are dropped.
- QC_INDIVIDUAL_MISSING_RATE (default: 0.1): Maximum missing rate per individual across genotypes. Individuals above the threshold are dropped.
- QC_HWE (default: 1e-50): Used to identify potential technical artifacts and drop variants.


### GWAS Parameters

- QC_MAF (default: 0.01): Minor Allele Frequency required to enter the GWAS pipeline.
- QC_MAC (default: 100): Minor Allele Count for REGENIE execution.


