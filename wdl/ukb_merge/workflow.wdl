version 1.0

struct BGENFileset {
    File bgen
    File bgi
    File sample
}

struct PLINKFileset {
    File bed
    File bim
    File fam
}

workflow merge_ukb_and_genomicc {
    input {
        String docker_image = "olivierlabayle/genomicc:main"
        Array[BGENFileset]+ bgen_filesets
        File hesin_critical_table
        File high_ld_regions
        PLINKFileset genomicc_genotypes
        PLINKFileset kgp_genotypes
        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"
        String qc_hwe_p = "1e-10"
        String qc_hwe_k = "0.001"
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
        String ancestry_threshold = "0.8"
        String palyndromic_threshold = "0.02"
        String julia_threads = "auto"
    }

    scatter (bgen_fileset in bgen_filesets) {
        call filter_ukb_chr {
            input:
                docker_image = docker_image,
                bgen_file = bgen_fileset.bgen,
                bgen_bgi_file = bgen_fileset.bgi,
                bgen_sample_file = bgen_fileset.sample,
                table_with_eids_to_exclude = hesin_critical_table,
                genomicc_genotyped_bim = genomicc_genotypes.bim,
                qc_genotype_missing_rate = qc_genotype_missing_rate,
                qc_individual_missing_rate = qc_individual_missing_rate
        }
    }

    scatter (set in filter_ukb_chr.plink_fileset) {
        File bed_files = set.bed
    }
    
    call merge_ukb_chrs {
        input:
            docker_image = docker_image,
            plink_filesets = filter_ukb_chr.plink_fileset,
            bed_files = bed_files
    }

    call align_ukb_variant_ids_with_kgp_and_keep_unrelated {
        input:
            docker_image = docker_image,
            ukb_bed = merge_ukb_chrs.merged_plink_fileset.bed,
            ukb_bim = merge_ukb_chrs.merged_plink_fileset.bim,
            ukb_fam = merge_ukb_chrs.merged_plink_fileset.fam,
            kgp_bed = kgp_genotypes.bed,
            kgp_bim = kgp_genotypes.bim,
            kgp_fam = kgp_genotypes.fam,
            palyndromic_threshold = palyndromic_threshold
    }

    call merge_genotypes_plink as merge_ukb_kgp {
        input:
            docker_image = docker_image,
            output_prefix = "ukb_kgp.merged",
            bed_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = kgp_genotypes.bed,
            bim_2 = kgp_genotypes.bim,
            fam_2 = kgp_genotypes.fam
    }

    call plink_qc as qc_ukb_kgp {
        input:
            docker_image = docker_image,
            bed_file = merge_ukb_kgp.merged_plink_fileset.bed,
            bim_file = merge_ukb_kgp.merged_plink_fileset.bim,
            fam_file = merge_ukb_kgp.merged_plink_fileset.fam,
            output_prefix = "ukb_kgp.merged.qc",
            qc_genotype_missing_rate = qc_genotype_missing_rate
    }

    call ld_prune as ld_prune_ukb_kgp {
        input:
            docker_image = docker_image,
            high_ld_regions = high_ld_regions,
            bed_file = qc_ukb_kgp.qc_fileset.bed,
            bim_file = qc_ukb_kgp.qc_fileset.bim,
            fam_file = qc_ukb_kgp.qc_fileset.fam,
            output_prefix = "ukb_kgp.merged.qc.ld_pruned",
            ip_values = ip_values,
            maf = maf
    }

    call estimate_ukb_ancestry_from_kgp {
        input:
            docker_image = docker_image,
            bed_file = ld_prune_ukb_kgp.ld_pruned_fileset.bed,
            bim_file = ld_prune_ukb_kgp.ld_pruned_fileset.bim,
            fam_file = ld_prune_ukb_kgp.ld_pruned_fileset.fam,
            output_filename = "ukb.ancestry_estimate.csv",
            ancestry_threshold = ancestry_threshold,
            cpus = julia_threads
    }

    call merge_genotypes_plink as merge_ukb_genomicc {
        input:
            docker_image = docker_image,
            output_prefix = "ukb_genomicc.merged",
            bed_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = genomicc_genotypes.bed,
            bim_2 = genomicc_genotypes.bim,
            fam_2 = genomicc_genotypes.fam
    }

    call plink_qc as qc_ukb_genomicc {
        input:
            docker_image = docker_image,
            bed_file = merge_ukb_genomicc.merged_plink_fileset.bed,
            bim_file = merge_ukb_genomicc.merged_plink_fileset.bim,
            fam_file = merge_ukb_genomicc.merged_plink_fileset.fam,
            output_prefix = "ukb_genomicc.merged.qc",
            qc_genotype_missing_rate = qc_genotype_missing_rate
    }

    output {
        PLINKFileset merged_ukb_fileset = merge_ukb_chrs.merged_plink_fileset
        PLINKFileset ukb_unrelated_fileset = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset
        PLINKFileset ukb_kgp_merged_fileset = merge_ukb_kgp.merged_plink_fileset
        PLINKFileset ukb_kgp_ld_pruned_fileset = ld_prune_ukb_kgp.ld_pruned_fileset
        File ancestry_estimate = estimate_ukb_ancestry_from_kgp.ancestry_estimate
        PLINKFileset ukb_genomicc_merged_fileset = merge_ukb_genomicc.merged_plink_fileset
        PLINKFileset ukb_genomicc_qc_fileset = qc_ukb_genomicc.qc_fileset
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

        plink2 \
            --max-alleles 2 \
            --pmerge-list merge_list.txt bfile \
            --output-chr chr26 \
            --make-bed \
            --out ukb_all_chr
    >>>

    output {
        PLINKFileset merged_plink_fileset = object {
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

task filter_ukb_chr {
    input {
        String docker_image
        File bgen_file
        File bgen_bgi_file
        File bgen_sample_file
        File table_with_eids_to_exclude
        File genomicc_genotyped_bim
        String qc_genotype_missing_rate
        String qc_individual_missing_rate
    }

    String bgen_prefix = basename(bgen_file, ".bgen")

    command <<<
        # Extract genomicc genotyped locations
        awk '{print $1, $4, $4}' ~{genomicc_genotyped_bim} > ranges_to_extract.txt

        # Extract sample IDs to exclude
        awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "eid") col=i; next} {print $col}' ~{table_with_eids_to_exclude} | sort -u > samples_to_remove.txt

        # Convert BGEN to PLINK, filtering both samples and variants
        plink2 \
            --bgen ~{bgen_file} ref-unknown --sample ~{bgen_sample_file} \
            --extract range ranges_to_extract.txt \
            --remove samples_to_remove.txt \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --max-alleles 2 \
            --make-bed \
            --out "~{bgen_prefix}.filtered"
    >>>

    output {
        PLINKFileset plink_fileset = object {
            bed: "${bgen_prefix}.filtered.bed",
            bim: "${bgen_prefix}.filtered.bim",
            fam: "${bgen_prefix}.filtered.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x8"
    }
}

task align_ukb_variant_ids_with_kgp_and_keep_unrelated {
    input {
        String docker_image
        File ukb_bed
        File ukb_bim
        File ukb_fam
        File kgp_bed
        File kgp_bim
        File kgp_fam
        String palyndromic_threshold = "0.02"
    }

    command <<<
        ukb_bed_prefix=$(dirname "~{ukb_bed}")/$(basename "~{ukb_bed}" .bed)
        kgp_bed_prefix=$(dirname "~{kgp_bed}")/$(basename "~{kgp_bed}" .bed)
        julia --project=/opt/sequential-gwas --startup-file=no --threads=auto /opt/sequential-gwas/bin/seq-gwas.jl \
            align-ukb-variant-ids-with-kgp-and-keep-unrelated \
            ${ukb_bed_prefix} \
            ${kgp_bed_prefix} \
            --threshold=~{palyndromic_threshold}
    >>>

    output {
        PLINKFileset ukb_unrelated_fileset = object {
            bed: "ukb_unrelated.bed",
            bim: "ukb_unrelated.bim",
            fam: "ukb_unrelated.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}

task merge_genotypes_plink {
    input {
        String docker_image
        String output_prefix = "merged_genotypes"
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
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset merged_plink_fileset = object {
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

task ld_prune {
    input {
        String docker_image
        File high_ld_regions
        File bed_file
        File bim_file
        File fam_file
        String output_prefix = "ld_pruned"
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
    }

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        plink2 \
            --bfile ${bed_prefix} \
            --indep-pairwise ~{ip_values}
        
        plink2 \
            --bfile ${bed_prefix} \
            --extract plink2.prune.in \
            --maf ~{maf} \
            --make-bed \
            --exclude range ~{high_ld_regions} \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x4"
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
        String cpus = "auto"
    }

    command <<<
        wget -O 20130606_g1k_3202_samples_ped_population.txt ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        julia --project=/opt/sequential-gwas --startup-file=no --threads=~{cpus} /opt/sequential-gwas/bin/seq-gwas.jl \
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
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}

task plink_qc {
    input {
        String docker_image
        File bed_file
        File bim_file
        File fam_file
        String output_prefix = "plink_qc"
        String qc_genotype_missing_rate = "0.02"
    }

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        plink2 \
            --bfile ${bed_prefix} \
            --geno ~{qc_genotype_missing_rate} \
            --make-bed \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset qc_fileset = object {
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}

# task export_critical_table { 
#     input {
#         File dataset
#         File fieldsfile
#     }

#     command {
#         dx extract_dataset ~{dataset} \
#         --fields-file ~{fieldsfile} \
#         -o="critical_table.csv" \
#         -icoding_option==RAW \
#         -iheader_style=UKB-FORMAT \
#         -ientity=hesin_critical \
#         -ifield_names_file_txt=~{fieldnames}
#     }

#     runtime {
#         dx_app: object {
#             id: "applet-xxxx",
#             type: "app" 
#         }
#         dx_timeout: "4H"
#         dx_instance_type: "mem1_ssd1_v2_x2"
#     }

#     output {
#         File outfile = "critical.csv"
#     }
# }