function cli_settings()
    s = ArgParseSettings(
        description="SeqGWAS",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(SequentialGWAS))
    )

    @add_arg_table! s begin
        "gwas-plots"
            action = :command
            help = "Generates GWAS plots."
        "pca-qc"
            action = :command
            help = "Runs PCAb-ased QC on genotypes to exclude outlier variants."

        "estimate-ancestry"
            action = :command
            help = "Estimate ancestry using the admixture software."

        "plot-pca"
            action = :command
            help = "Plots PCA of genotypes."

        "snps-to-flip"
            action = :command
            help = "Extracts SNPs to flip based on Illumina manifest file."

        "mock"
            action = :command
            help = "Creates mock data for testing purposes."

        "report-qc-effect"
            action = :command
            help = "Generates a report after a QC step."

        "qc-from-kgp"
            action = :command
            help = "Generates shared variants file and for each input, variants to flip, a new bim file and a summary of the suggested actions."
        
        "complete-bim-with-ref"
            action = :command
            help = "Completes a bim file with the KGP."
        
        "get-kgp-unrelated-individuals"
            action = :command
            help = "Completes a bim file with the KGP."

        "write-chromosomes"
            action = :command
            help = "Write list of chromosomes present in genotypes to file."

        "make-report"
            action = :command
            help = "Generates a report after the pipeline execution."
        
        "combine-covariates"
            action = :command
            help = "Merges covariates, ancestry and PCs files."

        "make-gwas-groups"
            action = :command
            help = "Generates groups, phenotypes and covariates files."
        
        "merge-covariates-pcs"
            action = :command
            help = "Merges covariates and PCs files."
    end

    @add_arg_table! s["gwas-plots"] begin
        "results"
            arg_type = String
            required = true
            help = "Path to GWAS results file."
        
        "group"
            arg_type = String
            required = true
            help = "Group name."
        
        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "gwas"
    end

    @add_arg_table! s["merge-covariates-pcs"] begin
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to covariates file."
        
        "pcs-file"
            arg_type = String
            required = true
            help = "Path to PCs file."

        "--output"
            arg_type = String
            help = "Output file name."
            default = "covariates_pcs.csv"
    end

    @add_arg_table! s["make-gwas-groups"] begin
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to covariates file."
        
        "variables-file"
            arg_type = String
            required = true
            help = "Path to variables file."

        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "group"
        
        "--min-group-size"
            arg_type = Int
            help = "Minimum group size."
            default = 100
    end

    @add_arg_table! s["combine-covariates"] begin
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to covariates file."
        
        "ancestry-file"
            arg_type = String
            required = true
            help = "Path to ancestry file."

        "pcs-file"
            arg_type = String
            required = true
            help = "Path to PCs file."

        "--output"
            arg_type = String
            help = "Output file name."
            default = "covariates.merged.csv"
    end

    @add_arg_table! s["make-report"] begin
        "unlifted_r8"
            arg_type = String
            required = true
            help = "Path to unlifted r8 file."

        "unlifted_2021_2023"
            arg_type = String
            required = true
            help = "Path to unlifted 2021_2023."

        "initial_bed_prefix_r8"
            arg_type = String
            required = true
            help = "Path to initial bed prefix r8."

        "initial_bed_prefix_2021_2023"
            arg_type = String
            required = true
            help = "Path to initial bed prefix 2021_2023."

        "initial_bed_prefix_2024_now"
            arg_type = String
            required = true
            help = "Path to initial bed prefix 2024_now."

        "basic_qc_prefix_r8"
            arg_type = String
            required = true
            help = "Path to basic qc prefix r8."

        "basic_qc_prefix_2021_2023"
            arg_type = String
            required = true
            help = "Path to basic qc prefix 2021_2023."

        "basic_qc_prefix_2024_now"
            arg_type = String
            required = true
            help = "Path to basic qc prefix 2024_now."

        "dup_samples_r8"
            arg_type = String
            required = true
            help = "Path to dup samples r8."

        "dup_samples_2021_2023"
            arg_type = String
            required = true
            help = "Path to dup samples 2021_2023."

        "dup-samples_2024_now"
            arg_type = String
            required = true
            help = "Path to dup samples 2024_now."

        "shared_variants"
            arg_type = String
            required = true
            help = "Path to shared variants file."

        "wgs_prefix"
            arg_type = String
            required = true
            help = "Path to WGS prefix."

        "merged_genotypes_prefix"
            arg_type = String
            required = true
            help = "Path to merged genotypes prefix."

        "unrelated_individuals"
            arg_type = String
            required = true
            help = "Path to unrelated individuals."

        "merged_qced_genotypes_prefix"
            arg_type = String
            required = true
            help = "Path to merged QCed genotypes prefix."

        "pca_plots_prefix"
            arg_type = String
            required = true
            help = "Path to PCA plots prefix."

        "high_loadings_variants"
            arg_type = String
            required = true
            help = "Path to high loadings variants."
            
        "final_genotypes_prefix"
            arg_type = String
            required = true
            help = "Prefix to final genotypes."
    end

    @add_arg_table! s["pca-qc"] begin
        "--input-prefix"
            arg_type = String
            required = true
            help = "Prefix to plink genotypes."
        "--output-prefix"
            arg_type = String
            required = true
            help = "Prefix to output files."
        "--ancestry-file"
            arg_type = String
            required = true
            help = "Path to ancestry file."
        "--npcs"
            arg_type = Int
            help = "Number of PCs to use."
            default = 10
        "--iqr-factor"
            arg_type = Float64
            help = "IQR factor to use for outlier detection."
            default = 3
        "--pca-approx"
            help = "Whether to use plink2 PCA approximation"
            action = :store_true
    end

    @add_arg_table! s["estimate-ancestry"] begin
        "genotypes-prefix"
            arg_type = String
            help = "Prefix to plink genotypes."
        
        "pedigree-file"
            arg_type = String
            help = "Path to pedigree file"

        "--output"
            arg_type = String
            help = "Output file name."
            default = "pop.pop"
    end

    @add_arg_table! s["plot-pca"] begin
        "eigenvectors"
            arg_type = String
            required = true
            help = "Path to eigenvectors file."
        
        "ancestry"
            arg_type = String
            required = true
            help = "Path to ancestry file."

        "--outprefix"
            arg_type = String
            help = "Prefix to output files."
            default = "pca"
    end

    @add_arg_table! s["complete-bim-with-ref"] begin
        "bim"
            arg_type = String
            help = "Path to bim file."

        "ref-bim"
            arg_type = String
            help = "Path to KGP bim file."
        
        "--output"
            arg_type = String
            help = "Output bim file."
            default = "complete.bim"
    end

    @add_arg_table! s["get-kgp-unrelated-individuals"] begin
        "pedigrees"
            arg_type = String
            required = true
            help = "Prefix to pedigree file."

        "--output"
            arg_type = String
            required = false
            default = "kgp_unrelated_individuals.txt"
            help = "Output file name."
    end

    @add_arg_table! s["write-chromosomes"] begin
        "input-prefix"
            arg_type = String
            required = true
            help = "Prefix to genotypes."

        "--output"
            arg_type = String
            required = false
            default = "chromosomes.txt"
            help = "Output file name."
    end

    @add_arg_table! s["qc-from-kgp"] begin
        "--outdir"
            arg_type = String
            required = false
            default = "."
            help = "Output directory."
        "--release-r8"
            arg_type = String
            required = true
            help = "Prefix corresponding to the release-r8."
        "--release-2021-2023"
            arg_type = String
            required = true
            help = "Prefix corresponding to the release-2021-2023."
        "--release-2024-now"
            arg_type = String
            required = true
            help = "Prefix corresponding to the release-2024-now."
        "--kgp"
            arg_type = String
            required = true
            help = "Prefix corresponding to the KGP."
        "--wgs-samples-file"
            arg_type = String
            required = true
            help = "Path to WGS sample IDs."
        "--threshold"
            arg_type = Float64
            help = "Palyndromic variants non-matching the 1000 GP alleles and whose MAF is close to 0.5 are dropped"
            default = 0.02
    end

    @add_arg_table! s["report-qc-effect"] begin
        "input-prefix"
            arg_type = String
            required = true
            help = "Prefix to input genotypes data."
        "output-prefix"
            arg_type = String
            required = true
            help = "Prefix to output genotypes data."
    end

    @add_arg_table! s["snps-to-flip"] begin
        "manifest-file"
            arg_type = String
            required = true
            help = "Path to Illumina manifest file."

        "--out"
            arg_type = String
            help = "Output .txt file with SNPs to flip (1 SNP per line)."
            default = "snps_to_flip.txt"
    end

    @add_arg_table! s["mock"] begin
        "release-r8"
            arg_type = String
            help = "Path prefix to the r8 release, i.e: `release-r8`(.ped | .map)."

        "release-2021-2023"
            arg_type = String
            help = "Path prefix to released data between 2021 and 2023, i.e: `release-2021-2023`(.ped | .map)."

        "release-after-2024"
            arg_type = String
            help = "Path prefix to released data after 2024, i.e: `release-after-2024`(.ped | .map)."

        "covariates"
            arg_type = String
            help = "Path to covariates file."
        
        "wgs-prefix"
            arg_type = String
            help = "Path prefix to WGS data."

        "--out-prefix"
            arg_type = String
            help = "Prefix to output data being generated."
            default = "mock"
        
        "--n-common-snps"
            arg_type = Int
            help = "Number of common SNPs to keep."
            default = 100
        
        "--n-distinct-snps"
            arg_type = Int
            help = "Number of distinct SNPs to keep."
            default = 10

        "--n-snps-wgs-not-genotyped"
            arg_type = Int
            help = "Number of SNPs in WGS data not genotyped in arrays."
            default = 30
        
        "--n-wgs-individuals"
            arg_type = Int
            help = "Number of WGS individuals."
            default = 10
            
        "--n-arrays-individuals"
            arg_type = Int
            help = "Number of arrays individuals."
            default = 1000
        
        "--rng"
            arg_type = Int
            help = "Random seed."
            default = 123

        "--verbosity"
            arg_type = Int
            help = "Verbosity level."
            default = 1
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    cmd_settings = settings[cmd]
    if cmd == "snps-to-flip"
        make_snps_to_flip_list(cmd_settings["out"], cmd_settings["manifest-file"])
    elseif cmd == "write-chromosomes"
        write_chromosomes(
            cmd_settings["input-prefix"], 
            output=cmd_settings["output"]
        )
    elseif cmd == "report-qc-effect"
        report_qc_effect(
            cmd_settings["input-prefix"], 
            cmd_settings["output-prefix"]
        )
    elseif cmd == "qc-from-kgp"
        generate_qc_extraction_files_from_kgp(
            cmd_settings["release-r8"],
            cmd_settings["release-2021-2023"],
            cmd_settings["release-2024-now"],
            cmd_settings["kgp"],
            cmd_settings["wgs-samples-file"];
            outdir=cmd_settings["outdir"],
            threshold=cmd_settings["threshold"]
        )
    elseif cmd == "complete-bim-with-ref"
        complete_bim_with_ref(
            cmd_settings["bim"],
            cmd_settings["ref-bim"],
            output=cmd_settings["output"]
        )
    elseif cmd == "get-kgp-unrelated-individuals"
        kgp_unrelated_individuals(
            cmd_settings["pedigrees"],
            outfile=cmd_settings["output"]
        )
    elseif cmd == "plot-pca"
        plot_pca(
            cmd_settings["eigenvectors"],
            cmd_settings["ancestry"];
            outprefix=cmd_settings["outprefix"]
        )
    elseif cmd == "estimate-ancestry"
        estimate_ancestry(
            cmd_settings["genotypes-prefix"],
            cmd_settings["pedigree-file"];
            output=cmd_settings["output"]
        )
    elseif cmd == "mock"
        mock_data(
            cmd_settings["release-r8"],
            cmd_settings["release-2021-2023"], 
            cmd_settings["release-after-2024"],
            cmd_settings["covariates"],
            cmd_settings["wgs-prefix"];
            outprefix=cmd_settings["out-prefix"],
            n_common_snps=cmd_settings["n-common-snps"], 
            n_distinct_snps=cmd_settings["n-distinct-snps"],
            n_snps_wgs_not_genotyped=cmd_settings["n-snps-wgs-not-genotyped"],
            n_wgs_individuals=cmd_settings["n-wgs-individuals"],
            n_arrays_individuals=cmd_settings["n-arrays-individuals"],
            rng=cmd_settings["rng"],
            verbosity=cmd_settings["verbosity"]
        )
    elseif cmd == "pca-qc"
        pca_qc(cmd_settings["input-prefix"],
            cmd_settings["ancestry-file"]; 
            npcs=cmd_settings["npcs"], 
            iqr_factor=cmd_settings["iqr-factor"],
            output_prefix = cmd_settings["output-prefix"],
            pca_approx=cmd_settings["pca-approx"]
        )
    elseif cmd == "make-report"
        make_report(;
            unlifted_r8 = cmd_settings["unlifted_r8"],
            unlifted_2021_2023 = cmd_settings["unlifted_2021_2023"],
            initial_bed_prefix_r8 = cmd_settings["initial_bed_prefix_r8"],
            initial_bed_prefix_2021_2023 = cmd_settings["initial_bed_prefix_2021_2023"],
            initial_bed_prefix_2024_now = cmd_settings["initial_bed_prefix_2024_now"],
            basic_qc_prefix_r8 = cmd_settings["basic_qc_prefix_r8"],
            basic_qc_prefix_2021_2023 = cmd_settings["basic_qc_prefix_2021_2023"],
            basic_qc_prefix_2024_now = cmd_settings["basic_qc_prefix_2024_now"],
            dup_samples_r8 = cmd_settings["dup_samples_r8"],
            dup_samples_2021_2023 = cmd_settings["dup_samples_2021_2023"],
            dup_samples_2024_now = cmd_settings["dup-samples_2024_now"],
            shared_variants = cmd_settings["shared_variants"],
            wgs_prefix = cmd_settings["wgs_prefix"],
            merged_genotypes_prefix = cmd_settings["merged_genotypes_prefix"],
            unrelated_individuals = cmd_settings["unrelated_individuals"],
            merged_qced_genotypes_prefix = cmd_settings["merged_qced_genotypes_prefix"],
            pca_plot_prefix = cmd_settings["pca_plots_prefix"],
            high_loadings_variants = cmd_settings["high_loadings_variants"],
            final_genotypes_prefix = cmd_settings["final_genotypes_prefix"]
        )
    elseif cmd == "combine-covariates"
        combine_covariates(
            cmd_settings["covariates-file"],
            cmd_settings["ancestry-file"],
            cmd_settings["pcs-file"];
            output=cmd_settings["output"]
        )
    elseif cmd == "make-gwas-groups"
        make_gwas_groups(
            cmd_settings["covariates-file"],
            cmd_settings["variables-file"];
            output_prefix=cmd_settings["output-prefix"],
            min_group_size=cmd_settings["min-group-size"]
        )
    elseif cmd == "merge-covariates-pcs"
        merge_covariates_pcs(
            cmd_settings["covariates-file"],
            cmd_settings["pcs-file"];
            output=cmd_settings["output"]
        )
    elseif cmd == "gwas-plots"
        gwas_plots(
            cmd_settings["results"],
            cmd_settings["group"];
            output_prefix=cmd_settings["output-prefix"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end