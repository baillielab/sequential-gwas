version 1.0

import "../common/structs.wdl"
import "../common/tasks.wdl"

struct RegenieStep1Files {
    Array[File] phenotypes_loco
    File list
}

workflow gwas {
    input {
        String docker_image = "olivierlabayle/genomicc:analysis_workflow"
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
        String approx_pca = "true"
        String maf = "0.01"
        String mac = "10"
        String ip_values = "1000 50 0.05"
        String regenie_cv_folds = "5"
        String regenie_bsize = "1000"
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

    # Merge covariates and LOCO PCs
    call merge_covariates_and_pcs {
        input:
            docker_image = docker_image,
            covariates_file = make_covariates_and_groups.updated_covariates,
            pcs_files = loco_pca.eigenvec,
            julia_cmd = get_julia_cmd.julia_cmd 
    }

    # Regenie Step 1
    scatter (group_individuals_and_plink_filesets in zip(make_covariates_and_groups.groups_lists, make_group_bed_qced.plink_fileset)) {

        File sample_list = group_individuals_and_plink_filesets.left
        PLINKFileset plink_fileset = group_individuals_and_plink_filesets.right

        call regenie_step1 {
            input:
                docker_image = docker_image,
                bed_file = plink_fileset.bed,
                bim_file = plink_fileset.bim,
                fam_file = plink_fileset.fam,
                sample_list = sample_list,
                covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                phenotypes_list = phenotypes,
                covariates_list = make_covariates_and_groups.covariates_list,
                cv_folds = regenie_cv_folds,
                bsize = regenie_bsize
        }
    }

    # Regenie Step 2

    scatter (pair in zip(make_covariates_and_groups.groups_lists, regenie_step1.step1_files)) {
        Pair[File, RegenieStep1Files] group_samples_and_regenie_step_1_files = pair
    }

    Array[Pair[PGENFileset, Pair[File, RegenieStep1Files]]] imputed_genotypes_to_groups_files = cross(imputed_genotypes, group_samples_and_regenie_step_1_files)

    scatter (imputed_genotype_and_group in imputed_genotypes_to_groups_files) {
        call regenie_step_2 {
            input:
                docker_image = docker_image,
                chr = imputed_genotype_and_group.left.chr,
                pgen_file = imputed_genotype_and_group.left.pgen,
                pvar_file = imputed_genotype_and_group.left.pvar,
                psam_file = imputed_genotype_and_group.left.psam,
                sample_list = imputed_genotype_and_group.right.left,
                covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                regenie_loco = imputed_genotype_and_group.right.right.phenotypes_loco,
                regenie_list = imputed_genotype_and_group.right.right.list,
                phenotypes_list = phenotypes,
                covariates_list = make_covariates_and_groups.covariates_list,
                npcs = npcs,
                bsize = regenie_bsize
        }
    }

    # Merge Regenie Step 2 results
    call merge_regenie_chr_results {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            regenie_step2_files = flatten(regenie_step_2.regenie_step2)
    }

    scatter (merged_results in merge_regenie_chr_results.merged_results) {
        # Generate GWAS plots
        call gwas_plots {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                results = merged_results
        }
    }
}

task gwas_plots {
    input {
        String docker_image
        String julia_cmd
        File results
    }

    command <<<
        ~{julia_cmd} gwas-plots \
            ~{results} \
            --output-prefix=gwas.plot
    >>>

    output {
        Array[File] plots = glob("gwas.plot*")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task merge_regenie_chr_results {
    input {
        String docker_image
        String julia_cmd
        Array[File] regenie_step2_files
    }

    command <<<
        for f in ~{sep=" " regenie_step2_files}; do
            echo "${f}"
        done > merge_list.txt

        ~{julia_cmd} merge-regenie-chr-results \
            merge_list.txt \
            --output-prefix=regenie.results
    >>>

    output {
        Array[File] merged_results = glob("regenie.results*.tsv")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task regenie_step_2 {
    input {
        String docker_image
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File sample_list
        File covariates_file
        Array[File] regenie_loco
        File regenie_list
        Array[String] phenotypes_list
        Array[String] covariates_list
        String npcs = "10"
        String bsize = "1000"
    }

    String group_name = sub(basename(sample_list, ".txt"), "grouped.individuals.", "")

    command <<<

        for file in  ~{sep=" " regenie_loco}; do
            ln -s "$file" .
        done

        input_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        pc_list=$(printf "CHR~{chr}_OUT_PC%s," {1..~{npcs}} | sed 's/,$//')
        full_covariates_list="~{sep="," covariates_list},${pc_list}"

        conda run -n regenie_env regenie \
            --step 2 \
            --pgen ${input_prefix} \
            --keep ~{sample_list} \
            --phenoFile ~{covariates_file} \
            --phenoColList ~{sep="," phenotypes_list} \
            --covarFile ~{covariates_file} \
            --covarColList ${full_covariates_list} \
            --bt \
            --firth --approx --pThresh 0.01 \
            --pred ~{regenie_list} \
            --bsize ~{bsize} \
            --out ~{group_name}.chr~{chr}.step2
    >>>

    output {
        Array[File] regenie_step2 = glob("${group_name}.chr${chr}.step2*.regenie")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}


task regenie_step1 {
    input {
        String docker_image
        File bed_file
        File bim_file
        File fam_file
        File sample_list
        File covariates_file
        Array[String] phenotypes_list
        Array[String] covariates_list
        String cv_folds = "5"
        String bsize = "1000"
    }

    String group_name = sub(basename(sample_list, ".txt"), "grouped.individuals.", "")

    command <<<
        genotypes_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        cv_option="--cv ~{cv_folds}"
        if [[ ~{cv_folds} == "loocv" ]]; then
            cv_option="--loocv"
        fi

        conda run -n regenie_env regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --keep ~{sample_list} \
            --phenoFile ~{covariates_file} \
            --phenoColList ~{sep="," phenotypes_list} \
            --covarFile ~{covariates_file} \
            --covarColList ~{sep="," covariates_list} \
            ${cv_option} \
            --bt \
            --bsize ~{bsize} \
            --lowmem \
            --out ~{group_name}.step1
        awk '{sub(".*/", "", $2); print $1, $2}' ~{group_name}.step1_pred.list > ~{group_name}.step1_pred.listrelative
    >>>

    output {
        RegenieStep1Files step1_files = object {
            phenotypes_loco: glob("${group_name}.step1_*.loco"),
            list: "${group_name}.step1_pred.listrelative"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

task merge_covariates_and_pcs {
    input {
        String docker_image
        File covariates_file
        Array[File] pcs_files
        String julia_cmd
    }

    command <<<
        for file in  ~{sep=" " pcs_files}; do
            ln -s "$file" .
        done

        ~{julia_cmd} \
            merge-covariates-pcs \
            ~{covariates_file} \
            pca \
            --output=merged_covariates_and_pcs.tsv
    >>>

    output {
        File covariates_and_pcs = "merged_covariates_and_pcs.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

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

    String output_prefix = "pca." + basename(bed_file, ".ldpruned.bed") + ".chr~{chr}_out"

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
        dx_instance_type: "mem2_ssd1_v2_x16"
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
        Array[String] covariates_list = read_lines("grouped.covariates_list.txt")
        Array[File] groups_lists = glob("grouped.individuals.*.txt")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
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
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}