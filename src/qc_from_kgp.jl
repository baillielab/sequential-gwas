
function kgp_minor_major_alleles_by_position_from_df(kgp_info)
    kgp_info_dict = Dict()
    for (variant_id, chr, pos, minor, major, minor_freq) in
            zip(kgp_info.VARIANT_ID, kgp_info.CHR_CODE, kgp_info.BP_COORD, kgp_info.MINOR_ALLELE, kgp_info.MAJOR_ALLELE, kgp_info.MINOR_ALLELE_FREQ)
        # Only SNPs are considered
        if length(minor) == length(major) == 1
            key = (chr, pos)
            # If there are multiple entries for the same position, these are multi-allelic variants that will be dropped
            if haskey(kgp_info_dict, key)
                kgp_info_dict[key] = nothing
            # Otherwise, we add the stats
            else
                kgp_info_dict[key] = (minor, major, minor_freq, variant_id)
            end
        end
    end
    return kgp_info_dict
end

function get_kgp_minor_major_alleles_by_position(kgp)
    kgp_info = GenomiccWorkflows.load_variants_info(kgp)
    return kgp_minor_major_alleles_by_position_from_df(kgp_info)
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

    kgp_minor, kgp_major, kgp_minor_freq, kgp_id = kgp_info[key]
    # Note: PLINK2 ALLELE_1 and ALLELE_2 are the minor and major alleles, respectively
    # Ideal case, where the variant in our dataset matches the KGP dataset
    if row.MINOR_ALLELE == kgp_minor && row.MAJOR_ALLELE == kgp_major
        return "KEEP (MINOR-MAJOR-MATCHING-KGP)"
    # There are two cases:
    ## (i) The variant is palindromic, then either:
    ##  - It wasn't properly flipped by GenomeStudio or earlier steps and needs to be flipped
    ##  - It has opposite allele frequency in our dataset (this should only happen when MAF ≈ 0.5 and we drop these cases)
    ## (ii) The variant is not palyndromic, it has an opposite allele frequency in our dataset (likely when MAF ≈ 0.5), we annotate it
    elseif row.MINOR_ALLELE == kgp_major && row.MAJOR_ALLELE == kgp_minor
        kgp_alleles = Set([kgp_minor, kgp_major])
        if kgp_alleles == Set(["A", "T"]) || kgp_alleles == Set(["C", "G"])
            if almost_balanced(row.MINOR_ALLELE_FREQ; threshold=threshold)
                return "DROP (PALINDROMIC-MINOR-MAJOR-REVERSED-KGP-BALANCED)"
            else
                return "FLIP (PALINDROMIC-MINOR-MAJOR-REVERSED-KGP-UNBALANCED)"
            end
        else
            return "KEEP (MINOR-MAJOR-REVERSED-KGP)" # Should be rare and mostly happen when MAF ≈ 0.5
        end
    else
        complement = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")
        # The variant alleles are the complement to the KGP alleles, we flip
        if row.MINOR_ALLELE == complement[kgp_minor] && row.MAJOR_ALLELE == complement[kgp_major]
            return "FLIP (COMPLEMENT-KGP)"
        # The variant alleles are the complement to the KGP alleles, but mismatch frequencies, we flip and annotate
        elseif row.MINOR_ALLELE == complement[kgp_major] && row.MAJOR_ALLELE == complement[kgp_minor]
            return "FLIP (COMPLEMENT-KGP-MINOR-MAJOR-REVERSED)"
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
            minor, major, minor_frequency, kgp_id = kgp_info[key]
            return kgp_id
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
        new_bim = DataFrames.select(bim, :CHR_CODE, :NEW_VARIANT_ID, :POSITION, :BP_COORD, :ALLELE_1, :ALLELE_2)
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
            DataFrames.select(to_flip, [:NEW_VARIANT_ID]),
            header=false
        )
    end
end

function update_bim_with_minor_major_info!(bim, acount)
        # Define a mapping from variant to minor and major alleles and minor allele frequency
    variant_id_to_minor_major = map(zip(acount.ID, acount.REF, acount.ALT, acount.ALT_CTS, acount.OBS_CT)) do (id, ref, alt, alt_ct, obs_ct)
        alt_freq = alt_ct / obs_ct
        minor, major, minor_freq = if alt_freq <= 0.5 
            (alt, ref, alt_freq) 
        else 
            (ref, alt, 1 - alt_freq)
        end
        id => (minor, major, minor_freq)
    end |> Dict
    # Update bim with MINOR_ALLELE, MAJOR_ALLELE and MINOR_ALLELE_FREQ
    bim.MINOR_ALLELE = [
        variant_id_to_minor_major[variant_id][1] for variant_id in bim.VARIANT_ID
    ]
    bim.MAJOR_ALLELE = [
        variant_id_to_minor_major[variant_id][2] for variant_id in bim.VARIANT_ID
    ]
    bim.MINOR_ALLELE_FREQ = [
        variant_id_to_minor_major[variant_id][3] for variant_id in bim.VARIANT_ID
    ]
    return bim
end

function load_variants_info(prefix)
    bim = read_bim(string(prefix, ".bim"))
    acount = CSV.read(string(prefix, ".acount"), DataFrame)
    return update_bim_with_minor_major_info!(bim, acount)
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
            DataFrames.select(matched_samples, [:FID, :IID]),
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
    release_r8_info = GenomiccWorkflows.load_variants_info(release_r8)
    release_2021_2023_info = GenomiccWorkflows.load_variants_info(release_2021_2023)
    release_2024_now_info = GenomiccWorkflows.load_variants_info(release_2024_now)
    kgp_info = GenomiccWorkflows.get_kgp_minor_major_alleles_by_position(kgp)

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
        release_2024_now = (basename(release_2024_now), read_fam(string(release_2024_now, ".fam"))),
        release_2021_2023 = (basename(release_2021_2023), read_fam(string(release_2021_2023, ".fam"))),
        release_r8 = (basename(release_r8), read_fam(string(release_r8, ".fam"))),
    )
    write_release_samples_to_drop(prefixes_and_fams, wgs_samples_file)
end