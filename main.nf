#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Input files params
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

include { CombineGenotypingArrays } from './workflows/combine_genotyping_arrays.nf'
include { ImputationWorkflow } from './workflows/imputation.nf'

workflow {
    genotyping_arrays = CombineGenotypingArrays()
}