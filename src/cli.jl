function cli_settings()
    s = ArgParseSettings(
        description="SeqGWAS",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(SequentialGWAS))
    )

    @add_arg_table! s begin
        "make-pop-file"
            action = :command
            help = "Creates a .pop file from a .fam and a pedigree file for the admixture software."

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
    end

    @add_arg_table! s["make-pop-file"] begin
        "fam-file"
            arg_type = String
            help = "Path to fam file"
        
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
        "--release-r8-bim"
            arg_type = String
            required = true
            help = "Bim file corresponding to the release-r8."
        "--release-2021-2023-bim"
            arg_type = String
            required = true
            help = "Bim file corresponding to the release-2021-2023."
        "--release-2024-now-bim"
            arg_type = String
            required = true
            help = "Bim file corresponding to the release-2024-now."
        "--kgp-bim"
            arg_type = String
            required = true
            help = "Bim file corresponding to the KGP."
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
        
        "--n-samples"
            arg_type = Int
            help = "Number of samples to keep."
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
            cmd_settings["release-r8-bim"],
            cmd_settings["release-2021-2023-bim"],
            cmd_settings["release-2024-now-bim"],
            cmd_settings["kgp-bim"],
            outdir=cmd_settings["outdir"],
            
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
            outprefix=cmd_settings["outprefix"]
        )
    elseif cmd == "make-pop-file"
        make_pop_file(
            cmd_settings["fam-file"],
            cmd_settings["pedigree-file"];
            output=cmd_settings["output"]
        )
    elseif cmd == "mock"
        mock_data(
            cmd_settings["release-r8"],
            cmd_settings["release-2021-2023"], 
            cmd_settings["release-after-2024"],
            cmd_settings["covariates"];
            outprefix=cmd_settings["out-prefix"],
            n_common_snps=cmd_settings["n-common-snps"], 
            n_distinct_snps=cmd_settings["n-distinct-snps"],
            n_samples=cmd_settings["n-samples"],
            rng=cmd_settings["rng"],
            verbosity=cmd_settings["verbosity"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end