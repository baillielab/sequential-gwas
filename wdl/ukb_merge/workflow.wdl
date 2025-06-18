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
        PLINKFileset genomicc_genotypes
        PLINKFileset kgp_genotypes
        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"
        String qc_hwe_p = "1e-10"
        String qc_hwe_k = "0.001"
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
                qc_individual_missing_rate = qc_individual_missing_rate,
                qc_hwe_p = qc_hwe_p,
                qc_hwe_k = qc_hwe_k
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
            palyndromic_threshold = "0.02"
    }

    output {
        PLINKFileset merged_ukb_fileset = merge_ukb_chrs.merged_plink_fileset
        PLINKFileset ukb_unrelated_fileset = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset
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
        String qc_hwe_p
        String qc_hwe_k
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
            --hwe ~{qc_hwe_p} ~{qc_hwe_k} \
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