version 1.0

import '../common/structs.wdl'

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
    String maf = "0.01"
    String mac = "10"
    }

    # Create groups and update covariates

    call make_covariates_and_groups {
        input:
            docker_image=docker_image,
            covariates_file=covariates_file,
            groupby=groupby,
            covariates=covariates,
            min_group_size=min_group_size,
            julia_use_sysimage=julia_use_sysimage,
            julia_threads=julia_threads
    }

    # Make BED files for each group

    scatter (sample_list in make_covariates_and_groups.groups_lists) {
        call make_group_bed_qced {
            input:
                docker_image = docker_image,
                genotypes_bed = genotypes.bed,
                genotypes_bim = genotypes.bim,
                genotypes_fam = genotypes.fam,
                sample_list = sample_list,
                maf = maf,
                mac = mac
        }
    }

  # scatter (i in range(length(make_group_bed_qced))) {
  #   call RegenieStep1 {
  #     input:
  #       phenotypes=make_covariates_and_groups.phenotypes,
  #       covariates=make_covariates_and_groups.covariates,
  #       covariates_list=make_covariates_and_groups.covariates_list,
  #       group="group${i}",
  #       samples=make_group_bed_qced[i].fam,
  #       bed=make_group_bed_qced[i].bed,
  #       bim=make_group_bed_qced[i].bim,
  #       fam=make_group_bed_qced[i].fam
  #   }
  # }

  # scatter (i in range(length(imputed_genotypes))) {
  #   call RegenieStep2 {
  #     input:
  #       group="group${i}",
  #       individuals=make_group_bed_qced[i].fam,
  #       step1_loco=RegenieStep1[i].loco,
  #       step1_pred=RegenieStep1[i].step1_pred,
  #       chr=chromosomes[i],
  #       imputed_genotypes=imputed_genotypes[i],
  #       covariates=make_covariates_and_groups.covariates,
  #       phenotypes=make_covariates_and_groups.phenotypes,
  #       covariates_list=make_covariates_and_groups.covariates_list
  #   }
  # }

  # call MergeRegenieResults {
  #   input:
  #     group="final",
  #     group_results=RegenieStep2.*.regenie
  # }

  # call MakeGWASPlots {
  #   input:
  #     group="final",
  #     group_results=MergeRegenieResults.output
  # }
}

task make_covariates_and_groups {
    input {
        String docker_image
        File covariates_file
        Array[String] groupby = []
        Array[String] covariates = ["AGE"]
        String min_group_size = "100"
        String julia_use_sysimage = "true"
        String julia_threads = "auto"
    }

    command <<<
        julia_cmd="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{julia_use_sysimage}" == "true" ]]; then
            julia_cmd+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{julia_threads}" == "auto" ]]; then
            julia_cmd+=" --threads=auto"
        fi

        groupby_string='~{sep="," groupby}'
        groupby_string_opt=""
        if [[ -n "${groupby_string}" ]]; then
            groupby_string_opt="--groupby=${groupby_string}"
        fi

        covariates_string='~{sep="," covariates}'

        ${julia_cmd} /opt/genomicc-workflows/bin/genomicc.jl \
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
        File bed = "${group_name}.bed"
        File bim = "${group_name}.bim"
        File fam = "${group_name}.fam"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

# task RegenieStep1 {
#   input {
#     File phenotypes
#     File covariates
#     File covariates_list
#     String group
#     File samples
#     File bed
#     File bim
#     File fam
#   }
#   command <<<
#     regenie \
#       --step 1 \
#       --bed ${bed.basename} \
#       --keep ${samples} \
#       --phenoFile ${phenotypes} \
#       --covarFile ${covariates} \
#       --covarColList $(cat ${covariates_list} | paste -sd, -) \
#       --cv 5 \
#       --bt \
#       --bsize 1000 \
#       --lowmem \
#       --out ${group}.step1
#     awk '{sub(".*/", "", $2); print $1, $2}' ${group}.step1_pred.list > ${group}.step1_pred.listrelative
#   >>>
#   output {
#     File loco = "${group}.step1_1.loco"
#     File step1_pred = "${group}.step1_pred.listrelative"
#   }
#   runtime {
#     cpu: 4
#     memory: "8G"
#   }
# }

# task RegenieStep2 {
#   input {
#     String group
#     File individuals
#     File step1_loco
#     File step1_pred
#     String chr
#     File imputed_genotypes
#     File covariates
#     File phenotypes
#     File covariates_list
#   }
#   command <<<
#     regenie \
#       --step 2 \
#       --pgen ${imputed_genotypes.basename} \
#       --phenoFile ${phenotypes} \
#       --covarFile ${covariates} \
#       --covarColList $(cat ${covariates_list} | paste -sd, -),CHR${chr}_OUT_PC1,CHR${chr}_OUT_PC2 \
#       --keep ${individuals} \
#       --bt \
#       --firth --approx --pThresh 0.01 \
#       --pred ${step1_pred} \
#       --bsize 1000 \
#       --out ${chr}.${group}.step2
#   >>>
#   output {
#     File regenie = "${chr}.${group}.step2.regenie"
#   }
#   runtime {
#     cpu: 4
#     memory: "16G"
#   }
# }

# task MergeRegenieResults {
#   input {
#     String group
#     Array[File] group_results
#   }
#   command <<<
#     ${get_julia_cmd()} merge-regenie-chr-results chr --output=${group}.csv
#   >>>
#   output {
#     File output = "${group}.csv"
#   }
# }

# task MakeGWASPlots {
#   input {
#     String group
#     File group_results
#   }
#   command <<<
#     ${get_julia_cmd()} gwas-plots ${group_results} ${group} --output-prefix ${group}
#   >>>
#   output {
#     File manhattan = "${group}.manhattan.png"
#     File qq = "${group}.qq.png"
#   }
# }