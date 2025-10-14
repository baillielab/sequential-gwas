#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Publish Directory
params.PUBLISH_DIR = "results"
params.KGP_PUBLISH_DIR = "${params.PUBLISH_DIR}/kgp"
params.ARRAY_GENOTYPES_PUBLISH_DIR = "${params.PUBLISH_DIR}/array_genotypes"
params.WGS_PUBLISH_DIR = "${params.PUBLISH_DIR}/wgs"
params.GATK_PUBLISH_DIR = "${params.PUBLISH_DIR}/gatk"
params.MERGED_PUBLISH_DIR = "${params.PUBLISH_DIR}/merged"
params.PCA_PUBLISH_DIR = "${params.PUBLISH_DIR}/pca"
params.ANCESTRY_PUBLISH_DIR = "${params.PUBLISH_DIR}/ancestry"

// Input files params
params.RESOURCES_DIR = "${projectDir}/assets/resources"
params.KGP_DIR = "${params.RESOURCES_DIR}/kgp"
params.GATK_DIR = "${params.RESOURCES_DIR}/gatk"
params.GRC37_TO_GRC38_CHAIN_FILE = "${projectDir}/assets/hg19ToHg38.over.chain"
params.HIGH_LD_REGIONS = "${projectDir}/assets/exclude_b38.txt" // Downloaded from https://github.com/GrindeLab/PCA/blob/main/data/highLD/exclude_b38.txt

// QC params
params.QC_GENOTYPE_MISSING_RATE = 0.02
params.QC_HWE_P = 1e-10
params.QC_HWE_K = 0.001
params.QC_INDIVIDUAL_MISSING_RATE = 0.02
params.IP_VALUES = "1000 50 0.05"
params.PCA_MAF = 0.01
params.PCA_APPROX = true
params.IQR_FACTOR = 10
params.FILTER_HIGH_LOADINGS_VARIANTS = false
params.ANCESTRY_THRESHOLD = 0.5
params.ANCESTRY_PROGRAM = "scope"
params.RELATEDNESS_DEGREE = 3

// Imputation
params.TOPMED_REFRESH_RATE = 180
params.TOPMED_MAX_PARALLEL_JOBS = 3
params.IMPUTATION_R2_FILTER = 0.9
params.N_SAMPLES_PER_IMPUTATION_JOBS = 19000
params.TOPMED_JOBS_LIST = "NO_TOPMED_JOBS"

// GWAS params
params.INFERRED_COVARIATES = "${projectDir}/assets/NO_INFERRED_COVARIATES"
params.MIN_GROUP_SIZE = 100
params.VARIABLES_CONFIG = "${projectDir}/assets/variables.yaml"
params.REGENIE_MAF = 0.01
params.REGENIE_MAC = 10
params.REGENIE_BSIZE = 1000
params.REGENIE_CV_NFOLDS = 5

// Other params
params.N_PCS = 20

include { CombineGeneticDatasets } from './workflows/combine_datasets.nf'
include { Imputation } from './workflows/imputation.nf'
include { KGP } from './workflows/kgp.nf'

log.info """\
         ${workflow.manifest.name} v${workflow.manifest.version}
         ==========================
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         profile      : ${workflow.profile}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}
         ==========================
         """
         .stripIndent()

workflow {
    CombineGeneticDatasets()
}