version 1.0

struct PGENFileset {
    String chr
    File pgen
    File psam
    File pvar
}

struct BGENFileset {
    String chr
    File bgen
    File bgi
    File sample
}

struct PLINKFileset {
    String chr
    File bed
    File bim
    File fam
}

struct BCFFileset {
    String chr
    File bcf
    File csi
}

workflow merge_ukb_and_genomicc {
    # Inputs

    input {
        String docker_image = "olivierlabayle/genomicc:main"
        Array[BGENFileset]+ bgen_filesets
        File hesin_critical_table
        File high_ld_regions
        PLINKFileset genomicc_genotypes
        Array[PGENFileset]+ genomicc_pgen_filesets
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
        String julia_use_sysimage = "true"
    }

    scatter (bgen_fileset in bgen_filesets) {
        call filter_ukb_chr {
            input:
                docker_image = docker_image,
                chr = bgen_fileset.chr,
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
            ukb_bed = merge_ukb_chrs.plink_fileset.bed,
            ukb_bim = merge_ukb_chrs.plink_fileset.bim,
            ukb_fam = merge_ukb_chrs.plink_fileset.fam,
            kgp_bed = kgp_genotypes.bed,
            kgp_bim = kgp_genotypes.bim,
            kgp_fam = kgp_genotypes.fam,
            palyndromic_threshold = palyndromic_threshold,
            julia_threads = julia_threads,
            julia_use_sysimage = julia_use_sysimage
    }

    call merge_genotypes_plink as merge_ukb_kgp {
        input:
            docker_image = docker_image,
            output_prefix = "ukb_kgp.merged",
            qc_genotype_missing_rate = qc_genotype_missing_rate,
            bed_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = kgp_genotypes.bed,
            bim_2 = kgp_genotypes.bim,
            fam_2 = kgp_genotypes.fam
    }

    call ld_prune as ld_prune_ukb_kgp {
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
            bed_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bed,
            bim_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.bim,
            fam_1 = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset.fam,
            bed_2 = genomicc_genotypes.bed,
            bim_2 = genomicc_genotypes.bim,
            fam_2 = genomicc_genotypes.fam
    }

    # Merging GenOMICC and UKB imputed genotypes

    scatter (bgen_fileset in bgen_filesets) {
        call bgen_to_vcf as ukb_bgen_to_vcf {
            input:
                docker_image = docker_image,
                output_prefix = "ukb.imputed",
                chr = bgen_fileset.chr,
                bgen_file = bgen_fileset.bgen,
                bgen_bgi_file = bgen_fileset.bgi,
                bgen_sample_file = bgen_fileset.sample,
                table_with_eids_to_exclude = hesin_critical_table,
                qc_genotype_missing_rate = qc_genotype_missing_rate,
                qc_individual_missing_rate = qc_individual_missing_rate
        }
    }

    scatter (pgen_fileset in genomicc_pgen_filesets) {
        call genomicc_pgen_to_bcf as genomicc_pgen_to_bcf {
            input:
                docker_image = docker_image,
                output_prefix = "genomicc.imputed",
                chr = pgen_fileset.chr,
                pgen_file = pgen_fileset.pgen,
                psam_file = pgen_fileset.psam,
                pvar_file = pgen_fileset.pvar
        }
    }

    scatter (ukb_genomicc_pair in zip(ukb_bgen_to_vcf.bcf_fileset, genomicc_pgen_to_bcf.bcf_fileset)) {
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

    # Outputs

    output {
        Array[PLINKFileset] filtered_ukb_chr = filter_ukb_chr.plink_fileset
        PLINKFileset merged_ukb_fileset = merge_ukb_chrs.plink_fileset
        PLINKFileset ukb_unrelated_fileset = align_ukb_variant_ids_with_kgp_and_keep_unrelated.ukb_unrelated_fileset
        PLINKFileset ukb_kgp_merged_fileset = merge_ukb_kgp.plink_fileset
        PLINKFileset ukb_kgp_ld_pruned_fileset = ld_prune_ukb_kgp.ld_pruned_fileset
        File ancestry_estimate = estimate_ukb_ancestry_from_kgp.ancestry_estimate
        PLINKFileset ukb_genomicc_merged_fileset = merge_ukb_genomicc.plink_fileset
        Array[BCFFileset] ukb_bcf_files = ukb_bgen_to_vcf.bcf_fileset
        Array[BCFFileset] genomicc_bcf_files = genomicc_pgen_to_bcf.bcf_fileset
        Array[PGENFileset] genomicc_ukb_pgen_filesets = merge_genomicc_ukb_bcfs_and_convert_to_pgen.pgen_fileset
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

task filter_ukb_chr {
    input {
        String docker_image
        String chr
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
        awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "eid") col=i; next} {print $col "\t" $col}' ~{table_with_eids_to_exclude} | sort -u > samples_to_remove.txt

        # Convert BGEN to PLINK, filtering both samples and variants
        plink2 \
            --bgen ~{bgen_file} ref-unknown --sample ~{bgen_sample_file} \
            --extract range ranges_to_extract.txt \
            --remove samples_to_remove.txt \
            --set-all-var-ids @:#:\$1:\$2 \
            --snps-only \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --output-chr chr26 \
            --max-alleles 2 \
            --make-bed \
            --out "~{bgen_prefix}.filtered"
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${bgen_prefix}.filtered.bed",
            bim: "${bgen_prefix}.filtered.bim",
            fam: "${bgen_prefix}.filtered.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x8"
    }

    meta {
        authors: ["Olivier Labayle"]
        description: "Filters a UKB BGEN chromosome file to be later merged with the GenOMICC genotyped dataset:\n- Variants: Only keep bi-allelic SNPs with genotype frequency higher than `qc_genotype_missing_rate` and matching the GenOMICC regions.\n- Samples: Exclude samples that have been critically ill and listed in `table_with_eids_to_exclude` and those with missing genotype rate higher than `qc_individual_missing_rate`."
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
        String julia_threads = "auto"
        String julia_use_sysimage = "true"
    }

    command <<<
        julia_cmd="julia --project=/opt/sequential-gwas --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/sequential-gwas/FlowOMMIC.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi
        ukb_bed_prefix=$(dirname "~{ukb_bed}")/$(basename "~{ukb_bed}" .bed)
        kgp_bed_prefix=$(dirname "~{kgp_bed}")/$(basename "~{kgp_bed}" .bed)
        ${julia_cmd} /opt/sequential-gwas/bin/seq-gwas.jl \
            align-ukb-variant-ids-with-kgp-and-keep-unrelated \
            ${ukb_bed_prefix} \
            ${kgp_bed_prefix} \
            --threshold=~{palyndromic_threshold}
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

task ld_prune {
    input {
        String docker_image
        File high_ld_regions
        String chr
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
            chr: chr,
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
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
        julia_cmd="julia --project=/opt/sequential-gwas --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/sequential-gwas/FlowOMMIC.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi

        wget -O 20130606_g1k_3202_samples_ped_population.txt ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        ${julia_cmd} /opt/sequential-gwas/bin/seq-gwas.jl \
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

task bgen_to_vcf {
    input {
        String docker_image
        String output_prefix
        String chr
        File bgen_file
        File bgen_bgi_file
        File bgen_sample_file
        File table_with_eids_to_exclude
        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"

    }

    command <<<
        # Extract sample IDs to exclude
        awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "eid") col=i; next} {print $col "\t" $col}' ~{table_with_eids_to_exclude} | sort -u > samples_to_remove.txt

        plink2 \
            --bgen ~{bgen_file} ref-unknown \
            --sample ~{bgen_sample_file} \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --output-chr chr26 \
            --export bcf \
            --out "~{output_prefix}.chr_~{chr}.temp"

        mamba run -n bcftools_env bcftools index "~{output_prefix}.chr_~{chr}.temp.bcf"
        mamba run -n bcftools_env bcftools norm \
            -m -both \
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

task genomicc_pgen_to_bcf {
    input {
        String docker_image
        String output_prefix
        String chr
        File pgen_file
        File psam_file
        File pvar_file
    }

    command <<<
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        plink2 \
            --pfile ${pgen_prefix} \
            --output-chr chr26 \
            --export bcf \
            --out "~{output_prefix}.chr_~{chr}.temp"
        
        mamba run -n bcftools_env bcftools index "~{output_prefix}.chr_~{chr}.temp.bcf"
        mamba run -n bcftools_env bcftools norm \
            -m -both \
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