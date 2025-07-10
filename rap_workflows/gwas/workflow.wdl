version 1.0

import '../common/structs.wdl'
import '../common/tasks.wdl'

workflow gwas {
    input {
        String docker_image = "olivierlabayle/genomicc:main"
        File covariates_file
        PLINKFileset genotypes
        Array[PGENFileset]+ imputed_genotypes
        Array[String] groupby = []
        Array[String] covariates = ["AGE", "SEX", "AGE_x_AGE", "AGE_x_SEX"]
        Array[String] phenotypes = ["SEVERE_COVID_19"]
        String min_group_size = "100"
        String julia_use_sysimage = "true"
        String julia_threads = "auto"
        File high_ld_regions = "assets/exclude_b38.txt"
        String npcs = "10"
        String approx_pca = "false"
        String maf = "0.01"
        String mac = "10"
        String ip_values = "1000 50 0.05"
    }

    call tasks.get_julia_cmd as get_julia_cmd {
        input:
            use_sysimage = julia_use_sysimage,
            threads = julia_threads
    }

    # Create groups and update covariates

    call make_covariates_and_groups {
        input:
            docker_image=docker_image,
            covariates_file=covariates_file,
            groupby=groupby,
            covariates=covariates,
            min_group_size=min_group_size,
            julia_cmd=get_julia_cmd.julia_cmd
    }

    # Make BED files for each group

    scatter (sample_list in make_covariates_and_groups.groups_lists) {
        call make_group_bed_qced {
            input:
                docker_image = docker_image,
                chr = genotypes.chr,
                genotypes_bed = genotypes.bed,
                genotypes_bim = genotypes.bim,
                genotypes_fam = genotypes.fam,
                sample_list = sample_list,
                maf = maf,
                mac = mac
        }
    }

    # LD prune BED files
    scatter (group_plink_fileset in make_group_bed_qced.plink_fileset) {
        call tasks.ld_prune as groups_ld_prune {
            input:
                docker_image = docker_image,
                high_ld_regions = high_ld_regions,
                chr = group_plink_fileset.chr,
                bed_file = group_plink_fileset.bed,
                bim_file = group_plink_fileset.bim,
                fam_file = group_plink_fileset.fam,
                output_prefix = basename(group_plink_fileset.bed, ".bed") + ".ldpruned",
                ip_values = ip_values,
                maf = maf
        }
    }

    # LOCO PCA on LD-pruned BED files
    scatter (imputed_genotype in imputed_genotypes) {
        String chromosomes = imputed_genotype.chr
    }

    Array[Pair[PLINKFileset, String]] groups_and_chrs = cross(groups_ld_prune.ld_pruned_fileset, chromosomes)

    scatter (group_and_chr in groups_and_chrs) {
        call loco_pca {
            input:
                docker_image = docker_image,
                chr = group_and_chr.right,
                bed_file = group_and_chr.left.bed,
                bim_file = group_and_chr.left.bim,
                fam_file = group_and_chr.left.fam,
                npcs = npcs,
                approx = approx_pca
            }
    }

}

# task merge_covariates_and_pcs {
#     input {
#         String docker_image
#         File covariates_file
#         Array[File] pcs_files
#     }

#     command <<<
#         ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
#             merge-covariates-and-pcs \
#             ~{covariates_file} \
#             ~{sep=" " pcs_files} \
#             --output-prefix=merged_covariates_and_pcs
#     >>>

#     output {
#         File merged_covariates = "merged_covariates_and_pcs.covariates.csv"
#         File pcs_list = "merged_covariates_and_pcs.pcs_list.txt"
#     }

#     runtime {
#         docker: docker_image
#         dx_instance_type: "mem1_ssd1_v2_x2"
#     }
# }

task loco_pca {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
        String npcs = 10
        String approx = "true"
    }

    String output_prefix = "pca" + basename(bed_file, ".ldpruned.bed") + ".chr~{chr}_out"

    command <<<
        genotypes_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        approx_option=""
        if [[ "~{approx}" == "true" ]]; then
            approx_option=" approx"
        fi

        plink2 \
            --bfile ${genotypes_prefix} \
            --not-chr ~{chr} \
            --pca ~{npcs}${approx_option} \
            --out ~{output_prefix}
    >>>

    output {
        File eigenvec = "${output_prefix}.eigenvec"
        File eigenval = "${output_prefix}.eigenval"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task make_covariates_and_groups {
    input {
        String docker_image
        File covariates_file
        Array[String] groupby = []
        Array[String] covariates = ["AGE"]
        String min_group_size = "100"
        String julia_cmd
    }

    command <<<
        groupby_string='~{sep="," groupby}'
        groupby_string_opt=""
        if [[ -n "${groupby_string}" ]]; then
            groupby_string_opt="--groupby=${groupby_string}"
        fi

        covariates_string='~{sep="," covariates}'

        ~{julia_cmd} \
            make-gwas-groups \
            ~{covariates_file} \
            --covariates=${covariates_string} \
            --output-prefix=grouped \
            --min-group-size=~{min_group_size} ${groupby_string_opt}
    >>>

    output {
        File updated_covariates = "grouped.covariates.csv"
        File covariates_list = "grouped.covariates_list.txt"
        Array[File] groups_lists = glob("grouped.individuals.*.txt")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task make_group_bed_qced {
    input {
        String docker_image
        String chr
        File genotypes_bed
        File genotypes_bim
        File genotypes_fam
        File sample_list
        String maf = "0.01"
        String mac = "10"
    }

    String group_name = sub(basename(sample_list, ".txt"), "grouped.individuals.", "")

    command <<<
        genotypes_prefix=$(dirname "~{genotypes_bed}")/$(basename "~{genotypes_bed}" .bed)

        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ~{sample_list} \
            --maf ~{maf} \
            --mac ~{mac} \
            --make-bed \
            --out ~{group_name}
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${group_name}.bed",
            bim: "${group_name}.bim",
            fam: "${group_name}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}