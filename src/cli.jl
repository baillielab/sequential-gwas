function cli_settings()
    s = ArgParseSettings(
        description="SeqGWAS",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(GenomiccWorkflows))
    )

    @add_arg_table! s begin
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
        
        "genotype-gvcf"
            action = :command
            help = "Completes a bim file with the KGP."
        
        "get-kgp-unrelated-individuals"
            action = :command
            help = "Completes a bim file with the KGP."

        "make-report"
            action = :command
            help = "Generates a report after the pipeline execution."
        
        "combine-covariates"
            action = :command
            help = "Merges covariates, ancestry and PCs files."
        
        "write-imputation-split-lists"
            action = :command
            help = "Writes imputation split lists for TOPMed API."
        
        "impute"
            action = :command
            help = "Imputes genotypes using the TOPMed API."

        "get-topmed-download-list"
            action = :command
            help = "Download results list using the TOPMed API."

        "download-topmed-file"
            action = :command
            help = "Downloads a file from TOPMed."
    end

    @add_arg_table! s["download-topmed-file"] begin
        "job-id"
            arg_type = String
            required = true
            help = "job-id"

        "token-file"
            arg_type = String
            required = true
            help = "Path to TOPMed API token file."

        "file-info"
            arg_type = String
            required = true
            help = "Info of the file to be downloaded."

        "--md5-file"
            arg_type = String
            default = nothing
            help = "Optional file to check MD5 against (for imputed .zip files)."

        "--refresh-rate"
            arg_type = Int
            help = "Rate at which to refresh the job status."
            default = 120
    end

    @add_arg_table! s["write-imputation-split-lists"] begin
        "genotypes-prefix"
            arg_type = String
            required = true
            help = "Prefix to genotypes"

        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "genomicc"

        "--n-samples-per-file"
            arg_type = Int
            help = "Number of samples per file."
            default = 20_000
    end

    @add_arg_table! s["get-topmed-download-list"] begin
        "job-id"
            arg_type = String
            required = true
            help = "job-id"

        "token-file"
            arg_type = String
            required = true
            help = "Path to TOPMed API token file."

        "--refresh-rate"
            arg_type = Int
            help = "Rate at which to refresh the job status."
            default = 120
    end

    @add_arg_table! s["impute"] begin
        "genotypes-prefix"
            arg_type = String
            required = true
            help = "Prefix to genotypes"

        "token-file"
            arg_type = String
            required = true
            help = "Path to TOPMed API token file."

        "--password"
            arg_type = String
            help = "Password for the TOPMed API."
            default = "abcde"

        "--max-concurrent-submissions"
            arg_type = Int
            help = "Maximum number of concurrent submissions to the TOPMed API."
            default = 3

        "--refresh-rate"
            arg_type = Int
            help = "Rate at which to refresh the job status."
            default = 120

        "--r2"
            arg_type = Float64
            help = "R2 threshold for imputation."
            default = 0.8

        "--samples-per-file"
            arg_type = Int
            help = "Number of samples per file."
            default = 10_000

        "--output-prefix"
            arg_type = String
            help = "Output prefix for jobs files."
            default = "."
    end

    @add_arg_table! s["combine-covariates"] begin
        "ancestry-file"
            arg_type = String
            required = true
            help = "Path to ancestry file."

        "pcs-file"
            arg_type = String
            required = true
            help = "Path to PCs file."
        
        "wgs-samples"
            arg_type = String
            required = true
            help = "Path to WGS samples file."
        
        "release-r8-fam"
            arg_type = String
            required = true
            help = "Path to r8 fam file."

        "release-2021-2023-fam"
            arg_type = String
            required = true
            help = "Path to release 2021-2023 fam file."
        
        "release-2024-now-fam"
            arg_type = String
            required = true
            help = "Path to release 2024-now fam file."

        "--output"
            arg_type = String
            help = "Output file name."
            default = "covariates.inferred.csv"
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

        "release_r8_qc_logs"
            arg_type = String
            required = true
            help = "Path to basic qc prefix r8."

        "release_2021_2023_qc_logs"
            arg_type = String
            required = true
            help = "Path to basic qc prefix 2021_2023."

        "release_2024_qc_logs"
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
        
        "release_r8_kgp_flip"
            arg_type = String
            required = true
            help = "Path to release r8 KGP flip file."

        "release_2021_2023_kgp_flip"
            arg_type = String
            required = true
            help = "Path to release 2021_2023 KGP flip file."

        "release_2024_now_kgp_flip"
            arg_type = String
            required = true
            help = "Path to release 2024_now KGP flip file."

        "shared_variants"
            arg_type = String
            required = true
            help = "Path to shared variants file."

        "wgs_prefix"
            arg_type = String
            required = true
            help = "Path to WGS prefix."

        "merge_log"
            arg_type = String
            required = true
            help = "Path to merged log file."

        "unrelated_individuals"
            arg_type = String
            required = true
            help = "Path to unrelated individuals."

        "qc_merge_log"
            arg_type = String
            required = true
            help = "Path to merged QCed log file."

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
        
        "covariates"
            arg_type = String
            required = true
            help = "Path to covariates file."
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

        "--threshold"
            arg_type = Float64
            help = "Threshold for ancestry assignment."
            default = 0.8

        "--program"
            arg_type = String
            help = "Ancestry estimation program (scope or admixture)."
            default = "scope"
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

    @add_arg_table! s["genotype-gvcf"] begin
        "gvcf-file"
            arg_type = String
            help = "Path to GVCF file."

        "shared-variants-plink"
            arg_type = String
            help = "Path to shared variants output by plink."

        "shared-variants-gatk"
            arg_type = String
            help = "Path to shared variants output for GATK (bed format)."

        "reference-genome"
            arg_type = String
            help = "Path to reference genome."
        
        "--output-prefix"
            arg_type = String
            help = "Output preix."
            default = "output"
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
    @info "Running GenOMICC Workflows: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "snps-to-flip"
        make_snps_to_flip_list(cmd_settings["out"], cmd_settings["manifest-file"])
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
    elseif cmd == "genotype-gvcf"
        genotype_gvcf(
            cmd_settings["gvcf-file"],
            cmd_settings["shared-variants-plink"],
            cmd_settings["shared-variants-gatk"],
            cmd_settings["reference-genome"];
            output_prefix=cmd_settings["output-prefix"]
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
            output=cmd_settings["output"],
            threshold=cmd_settings["threshold"],
            mode=cmd_settings["program"]
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
            release_r8_qc_logs = cmd_settings["release_r8_qc_logs"],
            release_2021_2023_qc_logs = cmd_settings["release_2021_2023_qc_logs"],
            release_2024_qc_logs = cmd_settings["release_2024_qc_logs"],
            dup_samples_r8 = cmd_settings["dup_samples_r8"],
            dup_samples_2021_2023 = cmd_settings["dup_samples_2021_2023"],
            dup_samples_2024_now = cmd_settings["dup-samples_2024_now"],
            release_r8_kgp_flip= cmd_settings["release_r8_kgp_flip"],
            release_2021_2023_kgp_flip = cmd_settings["release_2021_2023_kgp_flip"],
            release_2024_now_kgp_flip = cmd_settings["release_2024_now_kgp_flip"],
            shared_variants = cmd_settings["shared_variants"],
            wgs_prefix = cmd_settings["wgs_prefix"],
            merge_log = cmd_settings["merge_log"],
            unrelated_individuals = cmd_settings["unrelated_individuals"],
            qc_merge_log = cmd_settings["qc_merge_log"],
            pca_plot_prefix = cmd_settings["pca_plots_prefix"],
            high_loadings_variants = cmd_settings["high_loadings_variants"],
            final_genotypes_prefix = cmd_settings["final_genotypes_prefix"],
            covariates_file = cmd_settings["covariates"]
        )
    elseif cmd == "combine-covariates"
        combine_covariates(
            cmd_settings["ancestry-file"],
            cmd_settings["pcs-file"],
            cmd_settings["wgs-samples"],
            cmd_settings["release-r8-fam"],
            cmd_settings["release-2021-2023-fam"],
            cmd_settings["release-2024-now-fam"];
            output=cmd_settings["output"]
        )
    elseif cmd == "write-imputation-split-lists"
        write_imputation_split_lists(
            cmd_settings["genotypes-prefix"]; 
            output_prefix=cmd_settings["output-prefix"],
            samples_per_file=cmd_settings["n-samples-per-file"]
        )
    elseif cmd == "impute"
        impute(
            cmd_settings["genotypes-prefix"],
            cmd_settings["token-file"];
            password=cmd_settings["password"],
            max_concurrent_submissions=cmd_settings["max-concurrent-submissions"],
            refresh_rate=cmd_settings["refresh-rate"],
            r2=cmd_settings["r2"],
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "get-topmed-download-list"
        get_download_list_and_checksum(
            cmd_settings["job-id"], 
            cmd_settings["token-file"];
            refresh_rate=cmd_settings["refresh-rate"], 
        )
    elseif cmd == "download-topmed-file"
        download_topmed_file(
            cmd_settings["job-id"], 
            cmd_settings["token-file"],
            cmd_settings["file-info"];
            md5_file=cmd_settings["md5-file"],
            refresh_rate=cmd_settings["refresh-rate"]
            )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end