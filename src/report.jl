function get_filtered_variants_info(input_prefix, output_prefix)
    input_variants_info = read_bim(input_prefix)
    output_variants_info = read_bim(output_prefix)
    input_variants_info.VARIANT_ID = string.("chr", input_variants_info.CHR_CODE, ":", input_variants_info.BP_COORD)
    filtered_variants = DataFrame(VARIANT_ID = setdiff(input_variants_info.VARIANT_ID, output_variants_info.VARIANT_ID))
    return innerjoin(filtered_variants, input_variants_info, on=:VARIANT_ID)
end

function get_filtered_indiv_info(input_prefix, output_prefix)
    input_indiv_info = read_fam(input_prefix)
    output_indiv_info = read_fam(output_prefix)
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

function liftover_report(before_prefix, after_prefix)
    variants_info = DataFrame(
        read_map(string(before_prefix, ".map")), 
        [:CHR, :ID, :POS_BEFORE, :BP_BEFORE]
    )
    variants_after_liftover = DataFrame(
        read_map(string(after_prefix, ".map")), 
        [:CHR, :ID, :POS_AFTER, :BP_AFTER]
    )
    leftjoin!(variants_info, variants_after_liftover, on=[:CHR, :ID])
    variants_info.LIFTEDOVER = .!ismissing.(variants_info.BP_AFTER)
    CSV.write(string(before_prefix, ".liftover_report.csv"), variants_info)
end