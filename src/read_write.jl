
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