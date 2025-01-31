
"""
    read_bim(prefix)

Columns Description from: https://www.cog-genomics.org/plink/1.9/formats#bim
- Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
- Variant identifier
- Position in morgans or centimorgans (safe to use dummy value of '0')
- Base-pair coordinate (1-based; limited to 231-2)
- Allele 1 (corresponding to clear bits in .bed; usually minor)
- Allele 2 (corresponding to set bits in .bed; usually major)
"""
read_bim(prefix) = CSV.read(
    string(prefix, ".bim"), 
    DataFrame, 
    delim='\t', 
    header=["CHR_CODE", "VARIANT_ID", "POSITION", "BP_COORD", "ALLELE_1", "ALLELE_2"]
)

"""
    read_fam(prefix)

Columns Description from: https://www.cog-genomics.org/plink/1.9/formats#fam
- Family ID ('FID')
- Within-family ID ('IID'; cannot be '0')
- Within-family ID of father ('0' if father isn't in dataset)
- Within-family ID of mother ('0' if mother isn't in dataset)
- Sex code ('1' = male, '2' = female, '0' = unknown)
- Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
"""
read_fam(prefix) = CSV.read(
    string(prefix, ".fam"), 
    DataFrame, 
    delim=' ', 
    header=["FID", "IID", "FATHER_ID", "MOTHER", "SEX", "PHENOTYPE"]
)

function files_matching_prefix(prefix)
    directory, _prefix = splitdir(prefix)
    _directory = directory == "" ? "." : directory

    return map(
        f -> joinpath(directory, f),
        filter(
            f -> startswith(f, _prefix), 
            readdir(_directory)
        )
    )
end

function write_variants_intersection(output_file, input_dir)
    bim_files = filter(endswith(".bim"), readdir(input_dir))
    variants_intersection = mapreduce(intersect, bim_files) do file
        read_bim(joinpath(input_dir, file)[1:end-4]).VARIANT_ID
    end
    open(output_file, "w") do io
        for variant_id in variants_intersection
            println(io, variant_id)
        end
    end
end