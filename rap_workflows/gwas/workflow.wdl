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
        Array[String] filterby = []
        Array[String] covariates = ["AGE", "SEX", "AGE_x_AGE", "AGE_x_SEX"]
        Array[String] phenotypes = ["SEVERE_COVID_19"]
        String julia_use_sysimage = "true"
        String julia_threads = "auto"
        File high_ld_regions = "assets/exclude_b38.txt"
        # QC parameters
        String min_cases_controls = "10"
        String npcs = "10"
        String approx_pca = "true"
        String maf = "0.01"
        String mac = "10"
        String ip_values = "1000 50 0.05"
        # Regenie parameters
        String regenie_cv_folds = "5"
        String regenie_bsize = "1000"
        # Finemapping parameters
        String Xtype = "dosages" # or "genotypes"
        String min_sig_clump_size = "3"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "1e-5"
        String r2_threshold = "0.1"
        String clump_kb = "500"
        String n_causal = "10"
        String finemap_window_kb = "500"
        # Meta analysis
        Boolean meta_analysis = length(groupby) > 0
        Array[String] meta_exclude = ["ADMIXED"]
        String meta_method = "STDERR"
    }

    # Get generic Julia command
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
            filters=filterby,
            covariates=covariates,
            phenotypes_list=phenotypes,
            min_cases_controls=min_cases_controls,
            julia_cmd=get_julia_cmd.julia_cmd
    }


    # A GWAS is performed for each group defined by a combination of the user provided `groupby` and a phenotype
    scatter (sample_list in make_covariates_and_groups.groups_individuals) {
        String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

        # Make plink bed files for the group 
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

        # LD prune the bed file for LOCO PCA
        call tasks.ld_prune as groups_ld_prune {
            input:
                docker_image = docker_image,
                high_ld_regions = high_ld_regions,
                chr = make_group_bed_qced.plink_fileset.chr,
                bed_file = make_group_bed_qced.plink_fileset.bed,
                bim_file = make_group_bed_qced.plink_fileset.bim,
                fam_file = make_group_bed_qced.plink_fileset.fam,
                output_prefix = basename(make_group_bed_qced.plink_fileset.bed, ".bed") + ".ldpruned",
                ip_values = ip_values,
                maf = maf
        }

        # Perform LOCO PCA
        scatter (imputed_chr_fileset in imputed_genotypes) {
            call loco_pca {
                input:
                    docker_image = docker_image,
                    chr = imputed_chr_fileset.chr,
                    bed_file = groups_ld_prune.ld_pruned_fileset.bed,
                    bim_file = groups_ld_prune.ld_pruned_fileset.bim,
                    fam_file = groups_ld_prune.ld_pruned_fileset.fam,
                    npcs = npcs,
                    approx = approx_pca
            }
        }

        # Merge covariates and PCs
        call merge_covariates_and_pcs {
            input:
                docker_image = docker_image,
                group_name = group_name,
                covariates_file = make_covariates_and_groups.updated_covariates,
                pcs_files = loco_pca.eigenvec,
                julia_cmd = get_julia_cmd.julia_cmd 
        }

        # First run regenie step 1 using plink genotypes
        call regenie_step_1 {
            input:
                docker_image = docker_image,
                group_name = group_name,
                bed_file = make_group_bed_qced.plink_fileset.bed,
                bim_file = make_group_bed_qced.plink_fileset.bim,
                fam_file = make_group_bed_qced.plink_fileset.fam,
                sample_list = sample_list,
                covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                covariates_list = make_covariates_and_groups.covariates_list,
                cv_folds = regenie_cv_folds,
                bsize = regenie_bsize,
                maf = maf,
                mac = mac
        }

        # Second run regenie step 2 across imputed chromosomes filesets
        scatter (imputed_chr_fileset in imputed_genotypes) {
            # Perform Regenie Step 2
            call regenie_step_2 {
                input:
                    docker_image = docker_image,
                    group_name = group_name,
                    chr = imputed_chr_fileset.chr,
                    pgen_file = imputed_chr_fileset.pgen,
                    pvar_file = imputed_chr_fileset.pvar,
                    psam_file = imputed_chr_fileset.psam,
                    sample_list = sample_list,
                    covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                    regenie_loco = regenie_step_1.step1_files.phenotypes_loco,
                    regenie_list = regenie_step_1.step1_files.list,
                    covariates_list = make_covariates_and_groups.covariates_list,
                    npcs = npcs,
                    bsize = regenie_bsize,
                    mac = mac
            }

            # Finemap results for the chromosome
            call finemapping {
                input:
                    docker_image = docker_image,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    gwas_results = regenie_step_2.regenie_step2,
                    pgen_file = imputed_chr_fileset.pgen,
                    pvar_file = imputed_chr_fileset.pvar,
                    psam_file = imputed_chr_fileset.psam,
                    chr = imputed_chr_fileset.chr,
                    covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                    sample_file = regenie_step_2.regenie_step2_ids,
                    Xtype = Xtype,
                    group_name = group_name,
                    min_sig_clump_size = min_sig_clump_size,
                    lead_pvalue = lead_pvalue,
                    p2_pvalue = p2_pvalue,
                    r2_threshold = r2_threshold,
                    clump_kb = clump_kb,
                    n_causal = n_causal,
                    finemap_window_kb = finemap_window_kb
            }
        }

        # Merge GWAS results across chromosomes
        call merge_chr_results as merge_gwas_group_chr_results {
            input:
                docker_image = docker_image,
                output_prefix = group_name + ".gwas",
                julia_cmd = get_julia_cmd.julia_cmd,
                results_files = regenie_step_2.regenie_step2,
        }

        # Merge Finemapping results across chromosomes
        call merge_chr_results as merge_fp_group_chr_results {
            input:
                docker_image = docker_image,
                output_prefix = group_name + ".finemapping",
                julia_cmd = get_julia_cmd.julia_cmd,
                results_files = finemapping.finemapping_results
        }

        # Generate GWAS plots
        call gwas_plots as gwas_group_plots {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                gwas_results = merge_gwas_group_chr_results.merged_results,
                finemapping_results = merge_fp_group_chr_results.merged_results,
                maf = maf
        }
    }

    # Meta Analyse
    if (meta_analysis) {
        call meta_analyse {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                gwas_results = merge_gwas_group_chr_results.merged_results,
                exclude = meta_exclude,
                method = meta_method
        }

        # Finemap meta-analysed results: this effectively loops through phenotypes
        scatter (meta_gwas_result in meta_analyse.meta_gwas_results) {
            String phenotype = sub(basename(meta_gwas_result, ".gwas.tsv"), "META_ANALYSIS.", "")

            scatter (imputed_chr_fileset in imputed_genotypes) {
                call finemapping_summary_stats {
                    input:
                        docker_image = docker_image,
                        julia_cmd = get_julia_cmd.julia_cmd,
                        gwas_results = meta_gwas_result,
                        covariates_file = make_covariates_and_groups.updated_covariates,
                        sample_files = flatten(regenie_step_2.regenie_step2_ids),
                        pgen_file = imputed_chr_fileset.pgen,
                        pvar_file = imputed_chr_fileset.pvar,
                        psam_file = imputed_chr_fileset.psam,
                        chr = imputed_chr_fileset.chr,
                        min_sig_clump_size = min_sig_clump_size,
                        lead_pvalue = lead_pvalue,
                        p2_pvalue = p2_pvalue,
                        r2_threshold = r2_threshold,
                        clump_kb = clump_kb,
                        n_causal = n_causal,
                        finemap_window_kb = finemap_window_kb,
                        exclude = meta_exclude,
                        phenotype=phenotype
                }
            }

            call merge_chr_results as merge_fp_meta_chr_results {
                input:
                    docker_image = docker_image,
                    output_prefix = "META_ANALYSIS." + phenotype + ".finemapping" ,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    results_files = finemapping_summary_stats.finemapping_results
            }

            call gwas_plots as gwas_meta_plots {
                input:
                    docker_image = docker_image,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    gwas_results = meta_gwas_result,
                    finemapping_results = merge_fp_meta_chr_results.merged_results,
                    maf = maf
            }
        }

        Array[File] phenotypes_meta_plots = flatten(gwas_meta_plots.plots)
    }

    output {
        Array[File] gwas_group_results = merge_gwas_group_chr_results.merged_results
        Array[File] finemapping_group_results = merge_fp_group_chr_results.merged_results
        Array[File] group_plots = flatten(gwas_group_plots.plots)
        Array[File]? meta_gwas_results = meta_analyse.meta_gwas_results
        Array[File]? meta_finemapping_results = merge_fp_meta_chr_results.merged_results
        Array[File]? meta_plots = phenotypes_meta_plots
    }
}


task gwas_plots {
    input {
        String docker_image
        String julia_cmd
        File gwas_results
        File finemapping_results
        String maf = "0.01"
    }

    command <<<
        ~{julia_cmd} gwas-plots \
            ~{gwas_results} \
            ~{finemapping_results} \
            --maf=~{maf} \
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

task meta_analyse {
    input {
        String docker_image
        String julia_cmd
        Array[File] gwas_results
        Array[String] exclude
        String method = "STDERR"
    }

    command <<<
        for f in ~{sep=" " gwas_results}; do
            echo "${f}"
        done > gwas_meta_list.txt

        ~{julia_cmd} meta-analyse \
            gwas_meta_list.txt \
            --exclude=~{sep="," exclude} \
            --method=~{method} \
            --output-prefix=META_ANALYSIS
    >>>

    output {
        Array[File] meta_gwas_results = glob("META_ANALYSIS.*.tsv")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task merge_chr_results {
    input {
        String docker_image
        String output_prefix = "results.all_chr"
        String julia_cmd
        Array[File] results_files
    }

    command <<<
        for f in ~{sep=" " results_files}; do
            echo "${f}"
        done > merge_list.txt

        ~{julia_cmd} merge-chr-results \
            merge_list.txt \
            --output-prefix=~{output_prefix}
    >>>

    output {
        File merged_results = "${output_prefix}.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task finemapping_summary_stats {
    input {
        String docker_image
        String julia_cmd
        File gwas_results
        File pgen_file
        File pvar_file
        File psam_file
        String chr
        File covariates_file
        Array[File] sample_files
        String min_sig_clump_size = "3"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "1e-5"
        String r2_threshold = "0.1"
        String clump_kb = "1000"
        String n_causal = "10"
        String finemap_window_kb = "1000"
        Array[String] exclude = []
        String phenotype
    }

    command <<<
        # Make samples files list
        for f in ~{sep=" " sample_files}; do
            echo "${f}"
        done > samples_file.txt
        # Get PGEN prefix
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)
        # Run finemapping
        ~{julia_cmd} finemap \
            ~{gwas_results} \
            ${pgen_prefix} \
            ~{covariates_file} \
            samples_file.txt \
            --output-prefix=finemapping.meta_analysis.~{phenotype}.chr~{chr} \
            --min-sig-clump-size=~{min_sig_clump_size} \
            --lead-pvalue=~{lead_pvalue} \
            --p2-pvalue=~{p2_pvalue} \
            --r2-threshold=~{r2_threshold} \
            --clump-kb=~{clump_kb} \
            --n-causal=~{n_causal} \
            --finemap-window-kb=~{finemap_window_kb} \
            --phenotype=~{phenotype} \
            --rss \
            --exclude=~{sep="," exclude}
    >>>

    output {
        File finemapping_results = "finemapping.meta_analysis.${phenotype}.chr${chr}.tsv"
        File clumping_results = "finemapping.meta_analysis.${phenotype}.chr${chr}.clumps.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}


task finemapping {
    input {
        String docker_image
        String julia_cmd
        File gwas_results
        File pgen_file
        File pvar_file
        File psam_file
        String chr
        File covariates_file
        File sample_file
        String Xtype = "dosages" # or "genotypes"
        String group_name
        String min_sig_clump_size = "3"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "1e-5"
        String r2_threshold = "0.1"
        String clump_kb = "1000"
        String n_causal = "10"
        String finemap_window_kb = "1000"
    }

    command <<<
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        ~{julia_cmd} finemap \
            ~{gwas_results} \
            ${pgen_prefix} \
            ~{covariates_file} \
            ~{sample_file} \
            --Xtype=~{Xtype} \
            --output-prefix=finemapping.~{group_name}.chr~{chr} \
            --min-sig-clump-size=~{min_sig_clump_size} \
            --lead-pvalue=~{lead_pvalue} \
            --p2-pvalue=~{p2_pvalue} \
            --r2-threshold=~{r2_threshold} \
            --clump-kb=~{clump_kb} \
            --n-causal=~{n_causal} \
            --finemap-window-kb=~{finemap_window_kb} \
            --phenotype=${phenotype}
    >>>

    output {
        File finemapping_results = "finemapping.${group_name}.chr${chr}.tsv"
        File clumping_results = "finemapping.${group_name}.chr${chr}.clumps.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

task regenie_step_2 {
    input {
        String docker_image
        String group_name
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File sample_list
        File covariates_file
        Array[File] regenie_loco
        File regenie_list
        Array[String] covariates_list
        String npcs = "10"
        String bsize = "1000"
        String mac = "10"
    }

    command <<<

        for file in  ~{sep=" " regenie_loco}; do
            ln -s "$file" .
        done

        input_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        # Only retain bi-allelic (REGENIE can't handle non-biallelic)
        plink2 \
            --pfile ${input_prefix} \
            --keep ~{sample_list} \
            --max-alleles 2 \
            --min-alleles 2 \
            --rm-dup exclude-all list \
            --make-pgen \
            --out ${input_prefix}.biallelic_frequent
            
        # Make covariates list
        pc_list=$(printf "CHR~{chr}_OUT_PC%s," {1..~{npcs}} | sed 's/,$//')
        full_covariates_list="~{sep="," covariates_list},${pc_list}"

        # phenotype from group_name
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)
        echo $phenotype > phenotype.txt

        conda run -n regenie_env regenie \
            --step 2 \
            --pgen ${input_prefix}.biallelic_frequent \
            --keep ~{sample_list} \
            --phenoFile ~{covariates_file} \
            --phenoColList ${phenotype} \
            --write-samples \
            --covarFile ~{covariates_file} \
            --covarColList ${full_covariates_list} \
            --bt \
            --firth --approx --pThresh 0.01 \
            --minMAC ~{mac} \
            --pred ~{regenie_list} \
            --bsize ~{bsize} \
            --out ~{group_name}.chr~{chr}
    >>>

    output {
        File regenie_step2 = "${group_name}.chr${chr}_" + read_string("phenotype.txt") + ".regenie"
        File regenie_step2_ids = "${group_name}.chr${chr}_" + read_string("phenotype.txt") + ".regenie.ids"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}


task regenie_step_1 {
    input {
        String docker_image
        String group_name
        File bed_file
        File bim_file
        File fam_file
        File sample_list
        File covariates_file
        Array[String] covariates_list
        String cv_folds = "5"
        String bsize = "1000"
        String maf = "0.01"
        String mac = "10"
    }

    # String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

    command <<<
        genotypes_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        # In principle this is already taken care of in `make_group_bed_qced`, but we ensure here that the input is correct because it is cheap.
        # Only retain bi-allelic (REGENIE can't handle non-biallelic) and frequent variants within the sample list
        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ~{sample_list} \
            --mac ~{mac} \
            --maf ~{maf} \
            --min-alleles 2 \
            --max-alleles 2 \
            --write-snplist \
            --out biallelic_frequent

        # Parse cross-validation option
        cv_option="--cv ~{cv_folds}"
        if [[ ~{cv_folds} == "loocv" ]]; then
            cv_option="--loocv"
        fi

        # phenotype from group_name
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)

        conda run -n regenie_env regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --keep ~{sample_list} \
            --extract biallelic_frequent.snplist \
            --phenoFile ~{covariates_file} \
            --phenoColList ${phenotype} \
            --covarFile ~{covariates_file} \
            --covarColList ~{sep="," covariates_list} \
            --minMAC ~{mac} \
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
        String group_name
        File covariates_file
        Array[File] pcs_files
        String julia_cmd
    }

    String outfile = group_name + ".merged_covariates_and_pcs.tsv"

    command <<<
        for file in  ~{sep=" " pcs_files}; do
            ln -s "$file" .
        done

        ~{julia_cmd} \
            merge-covariates-pcs \
            ~{covariates_file} \
            pca \
            --output=~{outfile}
    >>>

    output {
        File covariates_and_pcs = "${outfile}"
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
        Array[String] filters = []
        Array[String] covariates = ["SEX", "AGE"]
        Array[String] phenotypes_list = ["SEVERE_COVID_19"]
        String min_cases_controls = "10"
        String julia_cmd
    }

    command <<<
        groupby_string='~{sep="," groupby}'
        groupby_string_opt=""
        if [[ -n "${groupby_string}" ]]; then
            groupby_string_opt="--groupby=${groupby_string}"
        fi

        filters_string='~{sep="," filters}'
        filters_string_opt=""
        if [[ -n "${filters_string}" ]]; then
            filters_string_opt="--filters=${filters_string}"
        fi

        covariates_string='~{sep="," covariates}'

        ~{julia_cmd} \
            make-gwas-groups \
            ~{covariates_file} \
            --covariates=${covariates_string} \
            --phenotypes=~{sep="," phenotypes_list} \
            --output-prefix=gwas \
            --min-cases-controls=~{min_cases_controls} ${groupby_string_opt} ${filters_string_opt}
    >>>

    output {
        File updated_covariates = "gwas.covariates.csv"
        Array[String] covariates_list = read_lines("gwas.covariates_list.txt")
        Array[File]+ groups_individuals = glob("gwas.individuals.*")
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

    String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

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