#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Publish Directory
params.PUBLISH_DIR = "results"
params.KGP_PUBLISH_DIR = "${params.PUBLISH_DIR}/kGP"
params.ARRAY_GENOTYPES_PUBLISH_DIR = "${params.PUBLISH_DIR}/array_genotypes"
params.WGS_PUBLISH_DIR = "${params.PUBLISH_DIR}/wgs"
params.GATK_PUBLISH_DIR = "${params.PUBLISH_DIR}/gatk"

// Input files params
params.RESOURCES_DIR = "${projectDir}/assets/resources/"
params.THOUSANDSGP_DIR = "${params.RESOURCES_DIR}/thousandsGP"
params.REFERENCE_GENOME = "${projectDir}/assets/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta" // From: https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
params.GRC37_TO_GRC38_CHAIN_FILE = "${projectDir}/assets/hg19ToHg38.over.chain.gz"
params.VARIANTS_TO_FLIP_GRC38 = "${projectDir}/assets/GSA-48v4-0_20085471_D2-minus-strand.txt"
params.VARIANTS_TO_FLIP_GRC37 = "${projectDir}/assets/GSA-24v3-0_A1-minus-strand.txt"

// QC params
params.QC_MAF = 0.01
params.QC_MAC = 100
params.QC_GENOTYPE_MISSING_RATE = 0.1
params.QC_HWE = 1e-50
params.QC_INDIVIDUAL_MISSING_RATE = 0.1

// Regenie params
params.REGENIE_BSIZE = 1000

include { AggregateGeneticData } from './workflows/aggregate_genetic_data.nf'
include { ImputationWorkflow } from './workflows/imputation.nf'
include { ThousandsGP } from './workflows/thousands_gp.nf'

workflow {
    AggregateGeneticData()
}