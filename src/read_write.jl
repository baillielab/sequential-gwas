
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
    read_fam(file)

Columns Description from: https://www.cog-genomics.org/plink/1.9/formats#fam
- Family ID ('FID')
- Within-family ID ('IID'; cannot be '0')
- Within-family ID of father ('0' if father isn't in dataset)
- Within-family ID of mother ('0' if mother isn't in dataset)
- Sex code ('1' = male, '2' = female, '0' = unknown)
- Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
"""
read_fam(file) = CSV.read(
    file, 
    DataFrame, 
    header=["FID", "IID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]
)

read_ped(file_prefix) = open(readlines, string(file_prefix, ".ped"))

write_ped(file_prefix, lines) = open(string(file_prefix, ".ped"), "w") do io 
    for line in lines
        println(io, line)
    end
end

"""
    read_map(file)

Columns Description from: https://www.cog-genomics.org/plink/1.9/formats

- Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
- Variant identifier
- Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
- Base-pair coordinate
"""
read_map(file) = CSV.read(file, DataFrame, delim='\t', header=["CHR", "ID", "POSITION", "BP_COORD"])

"""
    write_map(file_prefix, array)
"""
write_map(file_prefix, df) = CSV.write(string(file_prefix, ".map"), df, delim='\t', header=false)

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

function write_chromosomes(input_prefix; output="chromosomes.txt")
    chrs = SequentialGWAS.read_bim(string(input_prefix, ".bim")).CHR_CODE |> unique |> sort
    open(output, "w") do io
        for chr in chrs
            println(io, chr)
        end
    end
end

function complete_bim_with_ref(bim_file, ref_bim_file; output=bim_file)
    bim = SequentialGWAS.read_bim(bim_file)
    ref_bim = SequentialGWAS.read_bim(ref_bim_file)
    # Make the allele mapping for each variant
    variants_map = Dict()
    for row in Tables.namedtupleiterator(ref_bim)
        chr, pos, ref, alt = split(row.VARIANT_ID, ":")
        variants_map[string(chr, ":", pos)] = (
            variant_id = row.VARIANT_ID,
            allele_map = Dict(row.ALLELE_1 => row.ALLELE_2, row.ALLELE_2 => row.ALLELE_1)
        )
    end
    # Fill the missing allele and update the variant id
    bim.VARIANT_ID = convert(Vector{String}, bim.VARIANT_ID)
    bim.ALLELE_1 = convert(Vector{String}, bim.ALLELE_1)
    for row in eachrow(bim)
        info = variants_map[row.VARIANT_ID]
        row.VARIANT_ID = info.variant_id
        if row.ALLELE_1 == "."
            row.ALLELE_1 = info.allele_map[row.ALLELE_2]
        end
    end
    CSV.write(output, bim, header=false, delim="\t")
end

function make_pop_file(fam_file, pedigree_file; output="pop.pop")
    fam = SequentialGWAS.read_fam(fam_file)
    pedigree = CSV.read(pedigree_file, DataFrame, select=[:SampleID, :Superpopulation])
    leftjoin!(fam, pedigree, on=:IID => :SampleID)
    replace!(fam.Superpopulation, missing => "-")
    open(output, "w") do io
        for pop in fam.Superpopulation
            println(io, pop)
        end
    end
end