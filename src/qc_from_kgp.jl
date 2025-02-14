function set_action_column!(bim::DataFrame)
    bim.ACTION = map(eachrow(bim)) do row
        if row.VARIANT_ID === missing
            return "DROP (NOT-IN-KGP)"
        end
        if row.NON_BIALLELIC === true
            return "DROP (NON-BIALLELIC)"
        end
        # We now look at whether we need to flip the variant
        chr, pos, ref, alt = split(row.VARIANT_ID, ":")
        # Note: PLINK2 ALLELE_1 and ALLELE_2 are the minor and major alleles, respectively
        # Ideal case, where the variant in our dataset matches the KGP dataset
        if row.ALLELE_1 == alt && row.ALLELE_2 == ref
            return "KEEP (MATCHING-REF-ALT)"
        # There are two cases:
        ## (i) The variant is palyndromic, then either:
        ##  - It wasn't properly flipped by GenomeStudio (TODO: maybe GenomeStudio won't do that anymore)
        ##  - It has opposite allele frequency in our dataset (then all is fine)
        ## (ii) The variant is not palyndromic, it simply has an opposite allele frequency in our dataset (all is fine)
        ## Since opposite frequencies are not expected, we flip (i))
        elseif row.ALLELE_1 == ref && row.ALLELE_2 == alt
            alt_ref_set = Set([ref, alt])
            if alt_ref_set == Set(["A", "T"]) || alt_ref_set == Set(["C", "G"])
                return "FLIP (AMBIGUOUS-NOT-MATCHING-KGP)"
            else
                return "KEEP (REVERSE-REF-ALT)"
            end
        else
            complement = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")
            if (row.ALLELE_1 == complement[alt] && row.ALLELE_2 == complement[ref]) || 
                (row.ALLELE_1 == complement[ref] && row.ALLELE_2 == complement[alt])
                return "FLIP (COMPLEMENT)"
            else
                return "DROP (ALLELES-NOT-MATCHING-KGP)"
            end
        end
    end
end

"""
Compares the `bim` file to the `kgp_bim` and returns a new bim with an ACTION column
containing the potential following values:

1. DROP (REASON)
2. KEEP (REASON)
3. FLIP (REASON)

where "REASON" explains why the ACTION needs to be taken.
"""
function set_action_column(bim, kgp_bim)
    joined_bim = leftjoin(
        bim, 
        kgp_bim, 
        on=[:CHR_CODE, :BP_COORD]
    )
    select!(joined_bim, :CHR_CODE, 
        :KGP_VARIANT_ID => :VARIANT_ID, 
        :POSITION, 
        :BP_COORD, 
        :ALLELE_1, 
        :ALLELE_2
    )
    # Non unique variants are non-biallelic, we will drop them
    joined_bim.NON_BIALLELIC = nonunique(
        joined_bim, 
        [:CHR_CODE, :BP_COORD],
        keep = :noduplicates
    )
    joined_bim = unique(
        joined_bim, 
        [:CHR_CODE, :BP_COORD], 
        keep=:first
    )
    set_action_column!(joined_bim)
    return joined_bim
end

function generate_files_from_actions(outdir, prefixes_and_bims)
    # Get and Write variants present in all files.
    variants_intersection = intersect(
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_r8][2]).VARIANT_ID,
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_2021_2023][2]).VARIANT_ID,
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_2024_now][2]).VARIANT_ID
    )
    ## For plink
    CSV.write(
        joinpath(outdir, "variants_intersection.txt"), 
        DataFrame(VARIANT_ID=variants_intersection),
        header=false
    )
    ## For gatk
    gatk_shared_variants = DataFrame(
        mapreduce(x -> permutedims(x[1:2]), vcat, split.(variants_intersection, ":")),
        [:CHR_CODE, :BP_COORD]
    )
    gatk_shared_variants.BP_COORD = parse.(Int, gatk_shared_variants.BP_COORD) .- 1
    gatk_shared_variants.BP_COORD_END = gatk_shared_variants.BP_COORD .+ 1
    CSV.write(
        joinpath(outdir, "variants_intersection.bed"), 
        gatk_shared_variants, 
        delim='\t', 
        header=false
    )

    for (prefix, bim) in values(prefixes_and_bims)
        # Write summary
        CSV.write(
            joinpath(outdir, string(prefix, ".summary.csv")), 
            bim
        )
        # Write new bim file, replace missing by chr:bp:all1:all2, they will be dropped downstream
        new_bim = select(bim, :CHR_CODE, :VARIANT_ID, :POSITION, :BP_COORD, :ALLELE_1, :ALLELE_2)
        new_bim.VARIANT_ID = map(eachrow(new_bim)) do row
            if row.VARIANT_ID === missing
                return string(row.CHR_CODE, ":", row.BP_COORD, ":", row.ALLELE_1, ":", row.ALLELE_2)
            else
                return row.VARIANT_ID
            end
        end
        CSV.write(
            joinpath(outdir, string(prefix, ".new.bim")),
            new_bim, 
            header=false, 
            delim="\t"
        )
        # Write flip lists
        to_flip = filter(
            x -> startswith(x.ACTION, "FLIP") && x.VARIANT_ID in variants_intersection, 
            bim
        )
        CSV.write(
            joinpath(outdir, string(prefix, ".flip.txt")), 
            select(to_flip, [:VARIANT_ID]),
            header=false
        )
    end
end

function generate_qc_extraction_files_from_kgp(
    release_r8_bim_file, 
    release_2021_2023_bim_file, 
    release_2024_now_bim_file, 
    kgp_bim_file;
    outdir="."
    )
    # Load the KGP dataset
    release_r8_bim = SequentialGWAS.read_bim(release_r8_bim_file)
    release_2021_2023_bim = SequentialGWAS.read_bim(release_2021_2023_bim_file)
    release_2024_now_bim = SequentialGWAS.read_bim(release_2024_now_bim_file)
    kgp_bim = SequentialGWAS.read_bim(kgp_bim_file)
    rename!(kgp_bim, 
        :VARIANT_ID => :KGP_VARIANT_ID, 
        :POSITION => :KGP_POSITION,
        :ALLELE_1 => :KGP_ALLELE_1,
        :ALLELE_2 => :KGP_ALLELE_2,
        )
    # Format the chromosome and set the action column
    release_r8_bim = set_action_column(release_r8_bim, kgp_bim)
    release_2021_2023_bim = set_action_column(release_2021_2023_bim, kgp_bim)
    release_2024_now_bim = set_action_column(release_2024_now_bim, kgp_bim)
    # Write files to be used by plink
    prefixes_and_bims = Dict(
        :release_r8 => (first(splitext(basename(release_r8_bim_file))), release_r8_bim),
        :release_2021_2023 => (first(splitext(basename(release_2021_2023_bim_file))), release_2021_2023_bim),
        :release_2024_now => (first(splitext(basename(release_2024_now_bim_file))), release_2024_now_bim)
    )
    generate_files_from_actions(outdir, prefixes_and_bims)
end