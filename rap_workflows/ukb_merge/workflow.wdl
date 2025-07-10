version 1.0

import '../common/structs.wdl'
import '../common/tasks.wdl'

workflow merge_ukb_and_genomicc {
    # Inputs

    input {
        String docker_image = "olivierlabayle/genomicc:main"

        Array[BGENFileset]+ ukb_bgen_filesets
        File ukb_covariates
        File hesin_critical_table

        PLINKFileset genomicc_genotypes
        Array[PGENFileset]+ genomicc_pgen_filesets
        File genomicc_covariates
        File genomicc_inferred_covariates

        PLINKFileset kgp_genotypes

        File high_ld_regions
        File reference_genome
        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"
        String qc_hwe_p = "1e-10"
        String qc_hwe_k = "0.001"
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
        String ancestry_threshold = "0.8"
        String relatedness_degree = "3"
        String r2_threshold = "0.9"
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    # Index Reference Genome

    call index_reference_genome {
        input:
            docker_image = docker_image,
            reference_genome = reference_genome
    }

    # Filter UKB chromosomes with R2 and critical samples
    scatter (bgen_fileset in ukb_bgen_filesets) {
        call filter_ukb_chr_with_r2_and_critical_samples {
            input:
                docker_image = docker_image,
                chr = bgen_fileset.chr,
                bgen_file = bgen_fileset.bgen,
                bgen_bgi_file = bgen_fileset.bgi,
                bgen_sample_file = bgen_fileset.sample,
                vcf_info_file = bgen_fileset.vcf_info,
                vcf_info_file_index = bgen_fileset.vcf_info_index,
                table_with_eids_to_exclude = hesin_critical_table,
                qc_genotype_missing_rate = qc_genotype_missing_rate,
                qc_individual_missing_rate = qc_individual_missing_rate,
                r2_threshold = r2_threshold,
                julia_threads = julia_threads,
                julia_use_sysimage = julia_use_sysimage
        }
    }

    # Filter UKB BGEN files to keep only variants in GenOMICC

    scatter (pgen_fileset in filter_ukb_chr_with_r2_and_critical_samples.pgen_fileset) {
        call extract_genomicc_variants {
            input:
                docker_image = docker_image,
                chr = pgen_fileset.chr,
                pgen_file = pgen_fileset.pgen,
                pvar_file = pgen_fileset.pvar,
                psam_file = pgen_fileset.psam,
                genomicc_genotyped_bim = genomicc_genotypes.bim
        }
    }

    scatter (set in extract_genomicc_variants.plink_fileset) {
        File bed_files = set.bed
    }
    
    # Merge al lchromosomes together

    call merge_ukb_chrs {
        input:
            docker_image = docker_image,
            plink_filesets = extract_genomicc_variants.plink_fileset,
            bed_files = bed_files
    }

    # Make variant IDs consistent with KGP and keep only unrelated individuals

    call align_ukb_variants_with_kgp_and_keep_unrelated {
        input:
            docker_image = docker_image,
            ukb_bed = merge_ukb_chrs.plink_fileset.bed,
            ukb_bim = merge_ukb_chrs.plink_fileset.bim,
            ukb_fam = merge_ukb_chrs.plink_fileset.fam,
            kgp_bed = kgp_genotypes.bed,
            kgp_bim = kgp_genotypes.bim,
            kgp_fam = kgp_genotypes.fam,
            relatedness_degree = relatedness_degree,
            julia_threads = julia_threads,
            julia_use_sysimage = julia_use_sysimage
    }

    # Merging KGP genotypes with UKB genotypes for ancestry estimation

    call merge_genotypes_plink as merge_ukb_kgp {
        input:
            docker_image = docker_image,
            output_prefix = "ukb_kgp.merged",
            qc_genotype_missing_rate = qc_genotype_missing_rate,
            bed_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = kgp_genotypes.bed,
            bim_2 = kgp_genotypes.bim,
            fam_2 = kgp_genotypes.fam
    }

    # LD Pruning of the merged UKB and KGP genotypes

    call tasks.ld_prune as ld_prune_ukb_kgp {
        input:
            docker_image = docker_image,
            high_ld_regions = high_ld_regions,
            chr = merge_ukb_kgp.plink_fileset.chr,
            bed_file = merge_ukb_kgp.plink_fileset.bed,
            bim_file = merge_ukb_kgp.plink_fileset.bim,
            fam_file = merge_ukb_kgp.plink_fileset.fam,
            output_prefix = "ukb_kgp.merged.qc.ld_pruned",
            ip_values = ip_values,
            maf = maf
    }

    # Estimate ancestry from the LD-pruned UKB and KGP genotypes

    call estimate_ukb_ancestry_from_kgp {
        input:
            docker_image = docker_image,
            bed_file = ld_prune_ukb_kgp.ld_pruned_fileset.bed,
            bim_file = ld_prune_ukb_kgp.ld_pruned_fileset.bim,
            fam_file = ld_prune_ukb_kgp.ld_pruned_fileset.fam,
            output_filename = "ukb.ancestry_estimate.csv",
            ancestry_threshold = ancestry_threshold,
            julia_threads = julia_threads,
            julia_use_sysimage = julia_use_sysimage
    }

    # Merging GenOMICC genotypes with UKB matched genotypes (from imputed)

    call merge_genotypes_plink as merge_ukb_genomicc {
        input:
            docker_image = docker_image,
            output_prefix = "ukb_genomicc.merged",
            qc_genotype_missing_rate = qc_genotype_missing_rate,
            chr = "all",
            bed_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = genomicc_genotypes.bed,
            bim_2 = genomicc_genotypes.bim,
            fam_2 = genomicc_genotypes.fam
    }

    # Merging GenOMICC and UKB imputed genotypes

    scatter (pgen_fileset in filter_ukb_chr_with_r2_and_critical_samples.pgen_fileset) {
        call pgen_to_bcf as ukb_pgen_to_bcf {
            input:
                docker_image = docker_image,
                output_prefix = "ukb.imputed",
                chr = pgen_fileset.chr,
                pgen_file = pgen_fileset.pgen,
                pvar_file = pgen_fileset.pvar,
                psam_file = pgen_fileset.psam,
                reference_genome = reference_genome,
                reference_genome_index = index_reference_genome.reference_genome_index
        }
    }

    scatter (pgen_fileset in genomicc_pgen_filesets) {
        call pgen_to_bcf as genomicc_pgen_to_bcf {
            input:
                docker_image = docker_image,
                output_prefix = "genomicc.imputed",
                chr = pgen_fileset.chr,
                pgen_file = pgen_fileset.pgen,
                psam_file = pgen_fileset.psam,
                pvar_file = pgen_fileset.pvar,
                reference_genome = reference_genome,
                reference_genome_index = index_reference_genome.reference_genome_index
        }
    }

    scatter (ukb_genomicc_pair in zip(ukb_pgen_to_bcf.bcf_fileset, genomicc_pgen_to_bcf.bcf_fileset)) {
        call merge_genomicc_ukb_bcfs_and_convert_to_pgen {
            input:
                docker_image = docker_image,
                genomicc_chr = ukb_genomicc_pair.right.chr,
                genomicc_bcf = ukb_genomicc_pair.right.bcf,
                genomicc_csi = ukb_genomicc_pair.right.csi,
                ukb_chr = ukb_genomicc_pair.left.chr,
                ukb_bcf = ukb_genomicc_pair.left.bcf,
                ukb_csi = ukb_genomicc_pair.left.csi,
                qc_genotype_missing_rate = qc_genotype_missing_rate,
                output_prefix = "genomicc_ukb.merged.imputed"
        }
    }

    # Merge covariates

    call merge_genomicc_ukb_covariates {
        input:
            docker_image = docker_image,
            genomicc_covariates = genomicc_covariates,
            genomicc_inferred_covariates = genomicc_inferred_covariates,
            ukb_covariates = ukb_covariates,
            ukb_inferred_covariates = estimate_ukb_ancestry_from_kgp.ancestry_estimate,
            table_with_eids_to_exclude = hesin_critical_table,
            output_file = "ukb_genomicc.covariates.csv",
            julia_threads = julia_threads,
            julia_use_sysimage = julia_use_sysimage
    }

    # Make Report

    scatter (set in merge_genomicc_ukb_bcfs_and_convert_to_pgen.pgen_fileset) {
        File ukb_genomicc_merged_imputed_pvar_files = set.pvar
        File ukb_genomicc_merged_imputed_psam_files = set.psam
    }

    call make_report {
        input:
            docker_image = docker_image,
            ukb_kgp_merged_bed = merge_ukb_kgp.plink_fileset.bed,
            ukb_kgp_merged_bim = merge_ukb_kgp.plink_fileset.bim,
            ukb_kgp_merged_fam = merge_ukb_kgp.plink_fileset.fam,
            ukb_genomicc_merged_imputed_pvar = ukb_genomicc_merged_imputed_pvar_files,
            ukb_genomicc_merged_imputed_psam = ukb_genomicc_merged_imputed_psam_files,
            ukb_genomicc_merged_bim = merge_ukb_genomicc.plink_fileset.bim,
            ukb_genomicc_merged_fam = merge_ukb_genomicc.plink_fileset.fam,
            merged_covariates = merge_genomicc_ukb_covariates.merged_covariates,
            julia_threads = julia_threads,
            julia_use_sysimage = julia_use_sysimage
    }

    # Outputs

    output {
        Array[PGENFileset] filtered_ukb_r2_critical_pgen = filter_ukb_chr_with_r2_and_critical_samples.pgen_fileset
        Array[PLINKFileset] filtered_ukb_genomicc_variants_bed = extract_genomicc_variants.plink_fileset
        PLINKFileset merged_ukb_bed = merge_ukb_chrs.plink_fileset
        PLINKFileset ukb_unrelated_bed = align_ukb_variants_with_kgp_and_keep_unrelated.ukb_unrelated_fileset
        PLINKFileset ukb_kgp_merged_bed = merge_ukb_kgp.plink_fileset
        PLINKFileset ukb_kgp_ld_pruned_bed = ld_prune_ukb_kgp.ld_pruned_fileset
        File ukb_ancestry_estimates = estimate_ukb_ancestry_from_kgp.ancestry_estimate
        PLINKFileset ukb_genomicc_merged_bed = merge_ukb_genomicc.plink_fileset
        Array[BCFFileset] ukb_bcf_files = ukb_pgen_to_bcf.bcf_fileset
        Array[BCFFileset] genomicc_bcf_files = genomicc_pgen_to_bcf.bcf_fileset
        Array[PGENFileset] genomicc_ukb_merge_pgen = merge_genomicc_ukb_bcfs_and_convert_to_pgen.pgen_fileset
        File genomicc_ukb_merged_covariates = merge_genomicc_ukb_covariates.merged_covariates
        File report = make_report.report
    }
}

task merge_ukb_chrs {
    input {
        String docker_image
        Array[PLINKFileset] plink_filesets
        Array[File] bed_files
    }

    command <<<
        # Merge filesets
        for f in ~{sep=" " bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        plink \
            --biallelic-only \
            --merge-list merge_list.txt \
            --output-chr chr26 \
            --make-bed \
            --out ukb_all_chr
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: "all",
            bed: "ukb_all_chr.bed",
            bim: "ukb_all_chr.bim",
            fam: "ukb_all_chr.fam"
        }
        File merge_list = "merge_list.txt"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x4"
    }
}

task filter_ukb_chr_with_r2_and_critical_samples {
    input {
        String docker_image
        String chr
        File bgen_file
        File bgen_bgi_file
        File bgen_sample_file
        File vcf_info_file
        File vcf_info_file_index
        File table_with_eids_to_exclude
        String qc_genotype_missing_rate
        String qc_individual_missing_rate
        String r2_threshold = "0.9"
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    String input_prefix = basename(bgen_file, ".bgen")

    command <<<
        # Extract sample IDs to from table_with_eids_to_exclude
        awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "eid") col=i; next} {print $col "\t" $col}' ~{table_with_eids_to_exclude} | sort -u > samples_to_remove.txt
        # Convert to PGEN file
        plink2 \
            --bgen ~{bgen_file} ref-first \
            --sample ~{bgen_sample_file} \
            --output-chr chr26 \
            --make-pgen \
            --out ~{input_prefix}
        # Extract variant info from vcf file
        mamba run -n bcftools_env bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/R2\n' ~{vcf_info_file} > ~{input_prefix}.tsv
        # Check variants and get list of variants to extract
        julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi
        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
            make-ukb-bgen-qc-and-r2-filter-files \
            ~{input_prefix} \
            --threshold=~{r2_threshold} \
            --output=extract_list.txt
        # Filter PGEN file
        plink2 \
            --pfile ~{input_prefix} \
            --extract extract_list.txt \
            --remove samples_to_remove.txt \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --output-chr chr26 \
            --make-pgen \
            --out "~{input_prefix}.filtered"
    >>>

    output {
        PGENFileset pgen_fileset = object {
            chr: chr,
            pgen: "${input_prefix}.filtered.pgen",
            pvar: "${input_prefix}.filtered.pvar",
            psam: "${input_prefix}.filtered.psam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x16"
    }

}

task extract_genomicc_variants {
    input {
        String docker_image
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File genomicc_genotyped_bim
    }

    String pgen_prefix = basename(pgen_file, ".pgen")

    command <<<
        full_pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)
        # Extract genomicc genotyped locations
        awk '{print $1, $4, $4}' ~{genomicc_genotyped_bim} > ranges_to_extract.txt
        # Convert PGEN to PLINK, keep only bi-allelic SNPS, variant_IDS are reset to be unique to prevent multi-allelic variants obn multiple 
        # lines to cause problems during the merge. These are dropped later on and IDS set to match the 1000 GP ids.
        plink2 \
            --pfile ${full_pgen_prefix} \
            --extract range ranges_to_extract.txt \
            --set-all-var-ids @:#:\$1:\$2 \
            --snps-only \
            --output-chr chr26 \
            --max-alleles 2 \
            --make-bed \
            --out "~{pgen_prefix}.filtered"
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${pgen_prefix}.filtered.bed",
            bim: "${pgen_prefix}.filtered.bim",
            fam: "${pgen_prefix}.filtered.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x16"
    }

    meta {
        authors: ["Olivier Labayle"]
        description: "Filters a UKB BGEN chromosome file to be later merged with the GenOMICC genotyped dataset:\n- Variants: Only keep bi-allelic SNPs with genotype frequency higher than `qc_genotype_missing_rate` and matching the GenOMICC regions.\n- Samples: Exclude samples that have been critically ill and listed in `table_with_eids_to_exclude` and those with missing genotype rate higher than `qc_individual_missing_rate`."
    }
}

task align_ukb_variants_with_kgp_and_keep_unrelated {
    input {
        String docker_image
        File ukb_bed
        File ukb_bim
        File ukb_fam
        File kgp_bed
        File kgp_bim
        File kgp_fam
        String relatedness_degree = "3"
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    command <<<
        julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi
        ukb_bed_prefix=$(dirname "~{ukb_bed}")/$(basename "~{ukb_bed}" .bed)
        kgp_bed_prefix=$(dirname "~{kgp_bed}")/$(basename "~{kgp_bed}" .bed)
        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
            align-ukb-variants-with-kgp-and-keep-unrelated \
            ${ukb_bed_prefix} \
            ${kgp_bed_prefix} \
            --relatedness-degree=~{relatedness_degree}
    >>>

    output {
        PLINKFileset ukb_unrelated_fileset = object {
            chr: "all",
            bed: "ukb_unrelated.bed",
            bim: "ukb_unrelated.bim",
            fam: "ukb_unrelated.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
}

task merge_genotypes_plink {
    input {
        String docker_image
        String output_prefix = "merged_genotypes"
        String qc_genotype_missing_rate = "0.02"
        String chr = "all"
        File bed_1
        File bim_1
        File fam_1
        File bed_2
        File bim_2
        File fam_2
    }

    command <<<
        bed_prefix_1=$(dirname "~{bed_1}")/$(basename "~{bed_1}" .bed)
        bed_prefix_2=$(dirname "~{bed_2}")/$(basename "~{bed_2}" .bed)
        plink \
            --bfile ${bed_prefix_1} \
            --bmerge ${bed_prefix_2} \
            --output-chr chr26 \
            --biallelic-only strict \
            --make-bed \
            --out ~{output_prefix}.temp
        
        plink2 \
            --bfile ~{output_prefix}.temp \
            --geno ~{qc_genotype_missing_rate} \
            --make-bed \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x4"
    }
}

task estimate_ukb_ancestry_from_kgp {
    input {
        String docker_image
        File bed_file
        File bim_file
        File fam_file
        String output_filename = "ukb.ancestry_estimate.csv"
        String ancestry_threshold = "0.8"
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    command <<<
        julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi

        wget -O 20130606_g1k_3202_samples_ped_population.txt ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
            estimate-ancestry \
            ${bed_prefix} \
            20130606_g1k_3202_samples_ped_population.txt \
            --output=~{output_filename} \
            --threshold=~{ancestry_threshold}
    >>>

    output {
        File ancestry_estimate = output_filename
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x16"
    }
}

task pgen_to_bcf {
    input {
        String docker_image
        String output_prefix
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File reference_genome
        File reference_genome_index
    }

    command <<<
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        plink2 \
            --pfile ${pgen_prefix} \
            --fa ~{reference_genome} \
            --output-chr chr26 \
            --export bcf \
            --out "~{output_prefix}.chr_~{chr}.temp"

        mamba run -n bcftools_env bcftools index "~{output_prefix}.chr_~{chr}.temp.bcf"
        mamba run -n bcftools_env bcftools norm \
            -m -both \
            -f ~{reference_genome} \
            --check-ref wx \
            --output-type=b \
            --output="~{output_prefix}.chr_~{chr}.bcf" \
            --write-index=csi \
            "~{output_prefix}.chr_~{chr}.temp.bcf"
    >>>

    output {
        BCFFileset bcf_fileset = object {
            chr: chr,
            bcf: "${output_prefix}.chr_${chr}.bcf",
            csi: "${output_prefix}.chr_${chr}.bcf.csi"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x8"
    }
}

task merge_genomicc_ukb_bcfs_and_convert_to_pgen {
    input {
        String docker_image

        String genomicc_chr
        File genomicc_bcf
        File genomicc_csi

        String ukb_chr
        File ukb_bcf
        File ukb_csi

        String qc_genotype_missing_rate = "0.02"
        String output_prefix = "genomicc_ukb.merged.imputed"
    }

    command <<<
        if [[ "~{genomicc_chr}" != "~{ukb_chr}" ]]; then
            echo "GenOMICC and UKB chr files do not match. Make sure both input arrays have the same order." >&2
            exit 1
        fi

        # Merge BCF files
        mamba run -n bcftools_env bcftools merge \
            --output-type=b \
            --output="~{output_prefix}.chr_~{ukb_chr}.bcf" \
            --write-index=csi \
            ~{genomicc_bcf} ~{ukb_bcf}

        # Convert merged BCF to PGEN format
        plink2 \
            --bcf "~{output_prefix}.chr_~{ukb_chr}.bcf" \
            --geno ~{qc_genotype_missing_rate} \
            --make-pgen \
            --out "~{output_prefix}.chr_~{ukb_chr}"
    >>>

    output {
        PGENFileset pgen_fileset = object {
            chr: ukb_chr,
            pgen: "${output_prefix}.chr_~{ukb_chr}.pgen",
            psam: "${output_prefix}.chr_~{ukb_chr}.psam",
            pvar: "${output_prefix}.chr_~{ukb_chr}.pvar"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x8"
    }
}

task index_reference_genome {
    input {
        String docker_image
        File reference_genome
    }

    String reference_genome_filename = basename(reference_genome)

    command <<<
        mamba run -n bcftools_env samtools faidx ~{reference_genome}
        mv ~{reference_genome}.fai ~{reference_genome_filename}.fai
    >>>

    output {
        File reference_genome_index = "${reference_genome_filename}.fai"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}

task merge_genomicc_ukb_covariates {
    input {
        String docker_image
        File genomicc_covariates
        File genomicc_inferred_covariates
        File ukb_covariates
        File ukb_inferred_covariates
        File table_with_eids_to_exclude
        String output_file = "ukb_genomicc.covariates.csv"
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    command <<<
        julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi

        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
            merge-ukb-genomicc-covariates \
            ~{genomicc_covariates} \
            ~{genomicc_inferred_covariates} \
            ~{ukb_covariates} \
            ~{ukb_inferred_covariates} \
            ~{table_with_eids_to_exclude} \
            --output-file ~{output_file}
    >>>

    output {
        File merged_covariates = output_file
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}

task make_report {
    input {
        String docker_image

        File ukb_kgp_merged_bed
        File ukb_kgp_merged_bim
        File ukb_kgp_merged_fam

        Array[File] ukb_genomicc_merged_imputed_pvar
        Array[File] ukb_genomicc_merged_imputed_psam

        File ukb_genomicc_merged_bim
        File ukb_genomicc_merged_fam

        File merged_covariates

        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    command <<<
       julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi

        for f in ~{sep=" " ukb_genomicc_merged_imputed_pvar} ~{sep=" " ukb_genomicc_merged_imputed_psam}; do
            echo "${f%.bed}"
        done > ukb_genomicc_imputed_files_list.txt

        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
            make-ukb-genomicc-merge-report \
            ~{ukb_genomicc_merged_bim} \
            ~{ukb_genomicc_merged_fam} \
            ukb_genomicc_imputed_files_list.txt \
            ~{merged_covariates}
    >>>

    output {
        File report = "report.md"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}