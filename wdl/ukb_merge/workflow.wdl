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
        Array[BGENFileset]+ bgen_filesets
        File hesin_critical_table
        PLINKFileset genomicc_genotypes
        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"
        String qc_hwe_p = "1e-10"
        String qc_hwe_k = "0.001"
    }

    scatter (bgen_fileset in bgen_filesets) {
        call filter_ukb_chr {
            input: 
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
            plink_filesets = filter_ukb_chr.plink_fileset,
            bed_files = bed_files
    }

    output {
        PLINKFileset merged_ukb_fileset = merge_ukb_chrs.merged_plink_fileset
    }
}

task merge_ukb_chrs {
    input {
        Array[PLINKFileset] plink_filesets
        Array[File] bed_files
    }

    command <<<
        # Merge filesets
        for f in ~{sep=" " bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        plink2 \
            --pmerge-list merge_list.txt bfile \
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
        docker: "olivierlabayle/genomicc:main"
        dx_instance_type: "mem1_ssd2_v2_x4"
    }
}

task filter_ukb_chr {
    input {
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
        docker: "olivierlabayle/genomicc:main"
        dx_instance_type: "mem1_ssd2_v2_x8"
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