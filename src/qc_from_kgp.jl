
function get_kgp_ref_alt(kgp)
    kgp_info = SequentialGWAS.load_variants_info(kgp)
    kgp_ref_alt = Dict()
    for variant_id in kgp_info.VARIANT_ID
        chr, pos, ref, alt = split(variant_id, ":")
        if length(ref) == length(alt) == 1
            key = (chr, parse(Int, pos))
            if haskey(kgp_ref_alt, key)
                kgp_ref_alt[key] = nothing
            else
                kgp_ref_alt[key] = (ref, alt)
            end
        end
    end
    return kgp_ref_alt
end

almost_balanced(freq; threshold=0.02) = abs(freq - 0.5) < threshold

"""
Compares the variant to the 1000 GP information and returns an action to take together with a reason.

Potential actions are
1. DROP
2. KEEP
3. FLIP

Some particularly unexpected ACTION (REASON) are:

- "KEEP (REVERSE-REF-ALT)"
- "FLIP (COMPLEMENT-NOT-MATCHING-KGP)"

Because they mean the minor/major alleles are reversed in our dataset as compared to the reference KGP.
"""
function get_action(row, kgp_info; threshold=0.2)
    key = (row.CHR_CODE, row.BP_COORD)
    if !haskey(kgp_info, key)
        return "DROP (VARIANT-NOT-IN-KGP)"
    end
    
    if kgp_info[key] === nothing
        return "DROP (VARIANT-NON-BIALLELIC)"
    end

    ref, alt = kgp_info[key]
    # Note: PLINK2 ALLELE_1 and ALLELE_2 are the minor and major alleles, respectively
    # Ideal case, where the variant in our dataset matches the KGP dataset
    if row.ALLELE_1 == alt && row.ALLELE_2 == ref
        return "KEEP (MINOR-MAJOR-MATCHING-ALT-REF)"
    # There are two cases:
    ## (i) The variant is palindromic, then either:
    ##  - It wasn't properly flipped by GenomeStudio or earlier steps and needs to be flipped
    ##  - It has opposite allele frequency in our dataset (this should only happen when MAF ≈ 0.5 and we drop these cases)
    ## (ii) The variant is not palyndromic, it has an opposite allele frequency in our dataset (likely when MAF ≈ 0.5), we annotate it
    elseif row.ALLELE_1 == ref && row.ALLELE_2 == alt
        alt_ref_set = Set([ref, alt])
        if alt_ref_set == Set(["A", "T"]) || alt_ref_set == Set(["C", "G"])
            if almost_balanced(row.ALLELE_1_FREQ; threshold=threshold)
                return "DROP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF-BALANCED)"
            else
                return "FLIP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF)"
            end
        else
            return "KEEP (MINOR-MAJOR-REVERSED-ALT-REF)" # Should be rare and mostly happen when MAF ≈ 0.5
        end
    else
        complement = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")
        # The variant alleles are the complement to the KGP alleles, we flip
        if row.ALLELE_1 == complement[alt] && row.ALLELE_2 == complement[ref]
            return "FLIP (COMPLEMENT)"
        # The variant alleles are the complement to the KGP alleles, but mismatch frequencies, we flip and annotate
        elseif row.ALLELE_1 == complement[ref] && row.ALLELE_2 == complement[alt]
            return "FLIP (COMPLEMENT-MINOR-MAJOR-REVERSED-ALT-REF)"
        # All other cases are dropped (hopefully none)
        else
            return "DROP (ALLELES-NOT-MATCHING-KGP)"
        end
    end
end

function set_action_column!(variants_info::DataFrame, kgp_info; threshold=0.2)
    variants_info.ACTION = map(eachrow(variants_info)) do row
        get_action(row, kgp_info; threshold=threshold)
    end
end

function set_new_id_column!(variants_info::DataFrame, kgp_info)
    variants_info.NEW_VARIANT_ID = map(eachrow(variants_info)) do row
        key = (row.CHR_CODE, row.BP_COORD)
        if haskey(kgp_info, key) && kgp_info[key] !== nothing
            # Thew new id is made from chr:pos:ref:alt
            ref, alt = kgp_info[key]
            return string(row.CHR_CODE, ":", row.BP_COORD, ":", ref, ":", alt)
        else
            # The new id does not matter, this variant will be dropped
            return row.VARIANT_ID
        end
    end
end

function set_new_columns!(variants_info::DataFrame, kgp_info; threshold=0.2)
    set_action_column!(variants_info, kgp_info; threshold=threshold)
    set_new_id_column!(variants_info, kgp_info)
end

function generate_files_from_actions(outdir, prefixes_and_bims)
    # Get and write shared variants present in all files.
    variants_intersection = intersect(
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_r8][2]).NEW_VARIANT_ID,
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_2021_2023][2]).NEW_VARIANT_ID,
        filter(:ACTION => !startswith("DROP"), prefixes_and_bims[:release_2024_now][2]).NEW_VARIANT_ID
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
        new_bim = select(bim, :CHR_CODE, :NEW_VARIANT_ID, :POSITION, :BP_COORD, :ALLELE_1, :ALLELE_2)
        CSV.write(
            joinpath(outdir, string(prefix, ".new.bim")),
            new_bim, 
            header=false, 
            delim="\t"
        )
        # Write flip lists
        to_flip = filter(
            x -> startswith(x.ACTION, "FLIP") && x.NEW_VARIANT_ID in variants_intersection, 
            bim
        )
        CSV.write(
            joinpath(outdir, string(prefix, ".flip.txt")), 
            select(to_flip, [:NEW_VARIANT_ID]),
            header=false
        )
    end
end

function load_variants_info(prefix)
    bim = SequentialGWAS.read_bim(string(prefix, ".bim"))
    afreq = CSV.read(string(prefix, ".acount"), DataFrame)
    alt_freqs = Dict(
        zip(afreq.ID, afreq.ALT_CTS ./ afreq.OBS_CT)
    )
    bim.ALLELE_1_FREQ = map(bim.VARIANT_ID) do id
        get(alt_freqs, id, missing)
    end
    return bim
end

"""
We drop duplicate individuals according to the following priority:

WGS > More Recent Array > Older Array
"""
function write_release_samples_to_drop(prefixes_and_fams, wgs_samples_file)
    all_sample_ids = DataFrame(IID = readlines(wgs_samples_file))
    for (release_id, (release_prefix, fam_data)) in pairs(prefixes_and_fams)
        matched_samples = innerjoin(all_sample_ids, fam_data, on=:IID)
        CSV.write(
            string(release_prefix, ".samples_to_drop.txt"),
            select(matched_samples, [:FID, :IID]),
            header=false,
            delim="\t"
        )
        all_sample_ids = DataFrame(IID = union(all_sample_ids.IID, fam_data.IID))
    end
end

function generate_qc_extraction_files_from_kgp(
    release_r8, 
    release_2021_2023, 
    release_2024_now,
    kgp,
    wgs_samples_file;
    outdir=".",
    threshold=0.02
    )
    # Load the KGP dataset
    release_r8_info = SequentialGWAS.load_variants_info(release_r8)
    release_2021_2023_info = SequentialGWAS.load_variants_info(release_2021_2023)
    release_2024_now_info = SequentialGWAS.load_variants_info(release_2024_now)
    kgp_info = SequentialGWAS.get_kgp_ref_alt(kgp)

    # Format the chromosome and set the action column
    set_new_columns!(release_r8_info, kgp_info; threshold=threshold)
    set_new_columns!(release_2021_2023_info, kgp_info; threshold=threshold)
    set_new_columns!(release_2024_now_info, kgp_info; threshold=threshold)
    # Write variants related files to be used by plink:flip, drop
    prefixes_and_bims = Dict(
        :release_r8 => (basename(release_r8), release_r8_info),
        :release_2021_2023 => (basename(release_2021_2023), release_2021_2023_info),
        :release_2024_now => (basename(release_2024_now), release_2024_now_info)
    )
    generate_files_from_actions(outdir, prefixes_and_bims)
    # Write samples to drop from each release, the order is important here, it indicates the priority of the release
    prefixes_and_fams = (
        release_2024_now = (basename(release_2024_now), SequentialGWAS.read_fam(string(release_2024_now, ".fam"))),
        release_2021_2023 = (basename(release_2021_2023), SequentialGWAS.read_fam(string(release_2021_2023, ".fam"))),
        release_r8 = (basename(release_r8), SequentialGWAS.read_fam(string(release_r8, ".fam"))),
    )
    write_release_samples_to_drop(prefixes_and_fams, wgs_samples_file)
end