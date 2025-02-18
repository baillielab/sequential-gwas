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
params.GRC37_TO_GRC38_CHAIN_FILE = "${projectDir}/assets/hg19ToHg38.over.chain.gz"
params.VARIANTS_TO_FLIP_GRC38 = "${projectDir}/assets/GSA-48v4-0_20085471_D2-minus-strand.txt"
params.VARIANTS_TO_FLIP_GRC37 = "${projectDir}/assets/GSA-24v3-0_A1-minus-strand.txt"

// QC params
params.QC_GENOTYPE_MISSING_RATE = 0.1
params.QC_HWE = 1e-50
params.QC_INDIVIDUAL_MISSING_RATE = 0.1

// PCA Params
params.IP_VALUES = "1000 80 0.1"
params.PCA_MAF = 0.01
params.N_PCS = 10

// Regenie params
params.REGENIE_BSIZE = 1000

include { AggregateGeneticData } from './workflows/aggregate_genetic_data.nf'
include { ImputationWorkflow } from './workflows/imputation.nf'
include { KGP } from './workflows/kgp.nf'

workflow {
    AggregateGeneticData()
}