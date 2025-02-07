
"""
    read_bim(file)

Columns Description from: https://www.cog-genomics.org/plink/1.9/formats#bim
- Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
- Variant identifier
- Position in morgans or centimorgans (safe to use dummy value of '0')
- Base-pair coordinate (1-based; limited to 231-2)
- Allele 1 (corresponding to clear bits in .bed; usually minor)
- Allele 2 (corresponding to set bits in .bed; usually major)
"""
read_bim(file) = CSV.read(
    file, 
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
    header=["FID", "IID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]
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

function write_variants_intersection(output_prefix, input_dir)
    bim_files = filter(endswith(".bim"), readdir(input_dir))
    shared_variants = read_bim(joinpath(input_dir, popfirst!(bim_files)))
    # select!(shared_variants, [:CHR_CODE, :BP_COORD, :VARIANT_ID])
    for file in bim_files
        new_bim = read_bim(joinpath(input_dir, file))
        shared_variants = innerjoin(
            shared_variants, 
            select(new_bim, [:CHR_CODE, :BP_COORD]), 
            on=[:CHR_CODE, :BP_COORD]
        )
    end
    # Write bim
    CSV.write(
        string(output_prefix, ".bim"), 
        shared_variants, 
        delim='\t', 
        header=false
    )
    # Write plink file
    shared_variants.BP_COORD_END = shared_variants.BP_COORD
    CSV.write(
        string(output_prefix, ".csv"), 
        shared_variants[!, [:CHR_CODE, :BP_COORD, :BP_COORD_END, :VARIANT_ID]], 
        delim='\t', 
        header=false
    )
    # Write gatk file
    # shared_variants.BP_COORD_END .+= 1
    shared_variants.BP_COORD .-= 1
    CSV.write(
        string(output_prefix, ".bed"), 
        select(shared_variants, [:CHR_CODE, :BP_COORD, :BP_COORD_END]), 
        delim='\t', 
        header=false
    )
end

function write_chromosomes(input_prefix; output="chromosomes.txt")
    chrs = SequentialGWAS.read_bim(string(input_prefix, ".bim")).CHR_CODE |> unique |> sort
    open(output, "w") do io
        for chr in chrs
            println(io, chr)
        end
    end
end


function complete_bim_with_ref(bim_file, ref_bim_file)
    bim_file = SequentialGWAS.read_bim(bim_file)
    ref_bim_file = SequentialGWAS.read_bim(ref_bim_file)
    # Make the allele mapping for each variant
    variant_alleles_map = Dict{String, Dict{String, String}}()
    for row in Tables.namedtupleiterator(ref_bim_file)
        variant_alleles_map[row.VARIANT_ID] = Dict(row.ALLELE_1 => row.ALLELE_2, row.ALLELE_2 => row.ALLELE_1)
    end
    # Fill the missing allele
    for (idx, (variant_id, allele_1, allele_2)) in enumerate(zip(bim_file.VARIANT_ID, bim_file.ALLELE_1, bim_file.ALLELE_2))
        println(idx)
        if allele_1 == "."
            bim_file.ALLELE_1[idx] = variant_alleles_map[variant_id][allele_2]
        end
    end
end