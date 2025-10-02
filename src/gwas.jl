function merge_chr_results(merge_list_file; output_prefix = "results.all_chr")
    merge_list = readlines(merge_list_file)
    results = mapreduce(f -> CSV.read(f, DataFrame), vcat, merge_list)
    CSV.write(string(output_prefix, ".tsv"), results; delim="\t", header=true)
    return 0
end

get_chr_out_string(pc_filename) = splitext(splitext(pc_filename)[1])[2][2:end]

function read_loco_pcs(pc_file)
    chr_out = GenomiccWorkflows.get_chr_out_string(pc_file)
    pcs = CSV.read(pc_file, DataFrame, drop=["#FID"])
    PC_colnames = filter(!=("IID"), names(pcs))
    for PC_colname in PC_colnames
        rename!(pcs, Symbol(PC_colname) => Symbol(string(uppercase(chr_out), "_", PC_colname)))
    end
    return pcs
end

function merge_covariates_and_pcs(covariates_file, pcs_prefix; output="covariates_and_pcs.csv")
    covariates = CSV.read(covariates_file, DataFrame)
    pcs_dir = dirname(pcs_prefix)
    pcs_dir, pcs_prefix = pcs_dir == "" ? (".", "./$pcs_prefix") : (pcs_dir, pcs_prefix)
    pcs_files = filter(startswith(pcs_prefix), readdir(pcs_dir, join=true))
    chrs = unique(GenomiccWorkflows.get_chr_out_string.(pcs_files))
    chr_out_pcs = map(chrs) do chr
        chr_out_pcs_files = filter(x -> occursin(chr, x), pcs_files)
        mapreduce(read_loco_pcs, vcat, chr_out_pcs_files)
    end
    covariates_and_pcs = innerjoin(covariates, chr_out_pcs..., on=:IID)
    CSV.write(output, covariates_and_pcs, delim="\t", missingstring="NA")
    return 0
end

function make_gwas_groups(
    covariates_file; 
    groupby_string=nothing,
    covariates_string="AGE",
    phenotypes_string="SEVERE_COVID_19",
    filters_string=nothing,
    output_prefix="gwas", 
    min_cases_controls=100
    )
    phenotypes = split(phenotypes_string, ",")
    # Define additional covariates
    covariates, required_covariate_variables = read_and_process_covariates(covariates_file; covariates_string=covariates_string)
    # Write new covariates to file
    CSV.write(
        string(output_prefix, ".covariates.csv"), 
        covariates, 
        delim="\t"
    )
    # Write required covariates list to file for REGENIE
    open(string(output_prefix, ".covariates_list.txt"), "w") do io
        for covariate in required_covariate_variables
            println(io, covariate)
        end
    end
    # Make groups
    n_groups_passed = 0
    if groupby_string !== nothing
        groupby_variables = split(groupby_string, ",")
        for (groupkey, group) in pairs(groupby(covariates, groupby_variables, skipmissing=true, sort=true))
            group_id = join(groupkey, "_")
            n_phenotypes_passed = write_covariates_and_phenotypes_group(group, required_covariate_variables;
                group_id=group_id,
                phenotypes=phenotypes,
                output_prefix=output_prefix,
                min_cases_controls=min_cases_controls,
                filters_string=filters_string
            )
            n_groups_passed += n_phenotypes_passed
        end
    else
        n_groups_passed = write_covariates_and_phenotypes_group(covariates, required_covariate_variables; 
                group_id="all",
                phenotypes=phenotypes, 
                output_prefix=output_prefix, 
                min_cases_controls=min_cases_controls,
                filters_string=filters_string
        )
    end

    n_groups_passed > 0 || throw(ArgumentError("No group passed the min cases/controls threshold."))

    return 0
end

function group_and_phenotype_from_regenie_filename(filename)
    return split(replace(filename, 
        "results.all_chr." => "", 
        ".gwas.tsv" => ""
    ), ".")
end

group_needs_exclusion(group, exclude) =
    any(occursin(p, group) for p in exclude)

function run_metal_across_phenotypes!(regenie_files; output_prefix="gwas.meta_analysis", method="STDERR")
    tmp_dir = mktempdir()
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        phenotype = phenotype_key.PHENOTYPE
        metal_script = """
        # === DESCRIBE THE COLUMNS IN THE INPUT FILES ===
        MARKER ID 
        WEIGHT N 
        ALLELE ALLELE1 ALLELE0 
        FREQ A1FREQ 
        EFFECT BETA 
        STDERR SE 
        PVAL P_VAL
        SCHEME $method
        LOGPVALUE ON
        # === FOR EACH PHENOTYPE PROCESS / ANALYZE / ANALYZE HETEROGENEITY ===
        """
        for (group_file, group_basename) in zip(group.FILE, group.BASENAME)
            # Add P_VAL column expected by METAL
            group_gwas_results = CSV.read(group_file, DataFrame; delim="\t")
            transform!(group_gwas_results, :LOG10P => (x -> parse_pvalue.(x))  => :P_VAL)
            CSV.write(joinpath(tmp_dir, group_basename), group_gwas_results; delim="\t", header=true)
            # Add PROCESS command
            metal_script *= "PROCESS " * joinpath(tmp_dir, group_basename) * "\n"
        end
        metal_script *= "OUTFILE " * string(output_prefix, ".", phenotype, ". .tbl") * "\n"
        metal_script *= "ANALYZE HETEROGENEITY\n"
        group.METAL_FILE .= string(output_prefix, ".", phenotype, ".1.tbl")
        metal_script *= "QUIT"
        meta_script_file = joinpath(tmp_dir, "metal_script.$phenotype.txt")
        open(meta_script_file, "w") do io
            write(io, metal_script)
        end
        run(`metal $meta_script_file`)
    end
    return regenie_files
end

function post_process_metal_output(regenie_files; output_prefix="gwas.meta_analysis")
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        println(phenotype_key)
        metal_results = CSV.read(first(group.METAL_FILE), DataFrame; delim="\t")
        select!(metal_results, 
            "MarkerName" => "ID",
            "Effect" => "BETA",
            "StdErr" => "SE",
            "log(P)" => (x -> .-x) => "LOG10P", # -log10(P) is reported by Regenie, we make it compliant
            "Direction" => "DIRECTION",
            "HetISq" => "HET_ISQ",
            "HetChiSq" => "HET_CHISQ",
            "HetDf" => "HET_DF",
            "logHetP" => "LOG10P_HET"
        )
        append_GWAS_info_to_meta_analysis_results!(metal_results, group.FILE)
        CSV.write(string(output_prefix, ".", phenotype_key.PHENOTYPE, ".gwas.tsv"), metal_results; delim="\t", header=true)
    end
end

function append_GWAS_info_to_meta_analysis_results!(metal_results, phenotype_gwas_files)
    # Get variant info from GWAS results: CHROM, GENPOS, ALLELE0, ALLELE1, A1FREQ, N, NGROUPS
    variants_info_dict = Dict{String, Vector{Any}}()
    for phenotype_gwas_file in phenotype_gwas_files
        phenotype_gwas_results = CSV.read(phenotype_gwas_file, DataFrame; delim="\t")
        for row in Tables.namedtupleiterator(phenotype_gwas_results[!, [:CHROM, :GENPOS, :ID, :ALLELE0, :ALLELE1, :A1FREQ, :N]])
            if haskey(variants_info_dict, row.ID)
                variant_info = variants_info_dict[row.ID]
                variant_info[end] += 1 # update count of groups the variant was observed in
                variant_info[end-1] += row.N # update total N
                variant_info[end-2] = min(variant_info[end-2], row.A1FREQ) # update min A1FREQ
            else
                variants_info_dict[row.ID] = [row.CHROM, row.GENPOS, row.ALLELE0, row.ALLELE1, row.A1FREQ, row.N, 1]
            end
        end
    end
    # Update metal results with variant info
    transform!(metal_results, 
        :ID => ByRow(id -> variants_info_dict[id]) => [:CHROM, :GENPOS, :ALLELE0, :ALLELE1, :A1FREQ, :N, :NGROUPS]
    )

    return metal_results
end

function load_meta_analysis_worklist(regenie_files_list; exclude = [])
    regenie_files = DataFrame(FILE = readlines(regenie_files_list))
    regenie_files.BASENAME = basename.(regenie_files.FILE)
    transform!(regenie_files, 
        :BASENAME => ByRow(group_and_phenotype_from_regenie_filename) 
        => [:GROUP, :PHENOTYPE]
    )
    return filter!(:GROUP => (group -> !group_needs_exclusion(group, exclude)), regenie_files)
end

"""
    meta_analyse(regenie_files_list; exclude_string="ADMIXED", method="STDERR", output_prefix="gwas.meta_analysis")

Meta-analyse GWAS results from REGENIE using METAL. Groups and phenotypes are inferred from filenames. Results are meta analysed per phenotype across groups.

- regenie_files_list: path to a text file with a list of GWAS result files to meta-analyse
- exclude_string: comma-separated list of strings, if a group name contains any of these strings it will be excluded from meta-analysis (default: "ADMIXED")
- method: METAL meta-analysis method (default: "STDERR")
- output_prefix: prefix for output files. Per phenotype results are written to "<output_prefix>.<phenotype>.gwas.tsv". (default: "gwas.meta_analysis")
"""
function meta_analyse(regenie_files_list; exclude_string="ADMIXED", method="STDERR", output_prefix="gwas.meta_analysis")
    exclude = split(exclude_string, ",")
    regenie_files = load_meta_analysis_worklist(regenie_files_list; exclude = exclude)
    run_metal_across_phenotypes!(regenie_files; output_prefix=output_prefix, method=method)
    post_process_metal_output(regenie_files, output_prefix=output_prefix)

    return 0
end