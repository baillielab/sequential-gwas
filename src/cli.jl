function cli_settings()
    s = ArgParseSettings(
        description="SeqGWAS",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(SequentialGWAS))
    )

    @add_arg_table! s begin
        "snps-to-flip"
            action = :command
            help = "Extract SNPs to flip based on Illumina manifest file."

        "mock"
            action = :command
            help = "Create mock data for testing purposes."
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