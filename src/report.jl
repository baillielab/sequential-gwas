function get_filtered_variants_info(input_prefix, output_prefix)
    input_variants_info = read_bim(string(input_prefix, ".bim"))
    output_variants_info = read_bim(string(output_prefix, ".bim"))
    input_variants_info.VARIANT_ID = string.("chr", input_variants_info.CHR_CODE, ":", input_variants_info.BP_COORD)
    filtered_variants = DataFrame(VARIANT_ID = setdiff(input_variants_info.VARIANT_ID, output_variants_info.VARIANT_ID))
    return innerjoin(filtered_variants, input_variants_info, on=:VARIANT_ID)
end

function get_filtered_indiv_info(input_prefix, output_prefix)
    input_indiv_info = read_fam(string(input_prefix, ".fam"))
    output_indiv_info = read_fam(string(output_prefix, ".fam"))
    filtered_indiv = DataFrame(IID = setdiff(input_indiv_info.IID, output_indiv_info.IID))
    return innerjoin(filtered_indiv, input_indiv_info, on=:IID)
end

"""
    report_qc_effect(input_prefix, output_prefix)

Report the effect of the QC filtering process on SNPs and samples.
"""
function report_qc_effect(input_prefix, output_prefix)
    filtered_variants_info = get_filtered_variants_info(input_prefix, output_prefix)
    filtered_indiv_info = get_filtered_indiv_info(input_prefix, output_prefix)
    CSV.write(string(output_prefix, ".filtered_variants.csv"), filtered_variants_info)
    CSV.write(string(output_prefix, ".filtered_samples.csv"), filtered_indiv_info)
end



function make_report(;kwargs...)
    function prepend_args(script_string)
        args_string = join((string(key, " = \"", arg, "\" #hide") for (key, arg) in kwargs), "\n")
        return string(args_string, "\n", script_string)
    end
    Literate.markdown(
        joinpath(pkgdir(SequentialGWAS), "src", "report_template.jl"), 
        ".", 
        name="report", 
        flavor=Literate.CommonMarkFlavor(), 
        execute=true,
        preprocess=prepend_args
    )
end