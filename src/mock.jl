using CSV
using DataFrames
using DelimitedFiles
using Random

read_ped(file_prefix) = open(readlines, string(file_prefix, ".ped"))

write_ped(file_prefix, lines) = open(string(file_prefix, ".ped"), "w") do io 
    for line in lines
        println(io, line)
    end
end

read_map(file_prefix) = readdlm(string(file_prefix, ".map"))

write_map(file_prefix, array) = writedlm(string(file_prefix, ".map"), array, '\t')

outfile_prefix(prefix, file_id) = string(prefix, ".", file_id)

function make_mock_wgs(outprefix, wgs_prefix, wgs_snps, release_2024_now_map, new_sample_ids; 
    n_wgs_individuals=10,
    )
    # Make regions file to extract
    map_df = DataFrame(release_2024_now_map, ["CHR", "ID", "POS", "BP"])
    wgs_snps_df = DataFrame(ID=wgs_snps)
    merged = innerjoin(select(map_df, [:ID, :BP]), wgs_snps_df, on=:ID)
    tmpdir = mktempdir()
    regions_file = joinpath(tmpdir, "variants.tsv")
    CSV.write(variants_file, merged, delim='\t', header=false)
    # List gvcf files
    gvcf_files = filter(endswith("gvcf.gz"), readdir(dirname(wgs_prefix)))
    gvcf_files = rand(gvcf_files, n_wgs_individuals)
    # Extract regions from gvcf files
    for (input_gvcf, new_sample_id) in zip(gvcf_files, new_sample_ids.wgs)
        input_gvcf_filepath = joinpath(wgs_prefix, input_gvcf)
        output_gvcf_filepath = string(outprefix, ".", basename(input_gvcf))
        tmpdir = mktempdir()
        new_sample_file = joinpath(tmpdir, "new_sample.txt")
        write(new_sample_file, new_sample_id)
        run(`bcftools view -O z -R $regions_file $input_gvcf_filepath \| bcftools reheader -s $new_sample_file -o $output_gvcf_filepath`)
    end
end

function make_mock_map_files_and_wgs(
    release_r8, 
    release_2021_2023, 
    release_2024_now, 
    wgs_prefix,
    new_sample_ids; 
    outprefix="mock",
    n_common_snps=100,
    n_distinct_snps=10,
    n_wgs_individuals=10,
    n_snps_wgs_not_genotyped=30
    )
    origin_map_files = (
        release_r8 = read_map(release_r8),
        release_2021_2023 = read_map(release_2021_2023),
        release_2024_now = read_map(release_2024_now)
    )

    snp_intersection = intersect((mapfile[:, 2] for mapfile in origin_map_files)...)
    wgs_snps= rand(snp_intersection, n_common_snps+n_snps_wgs_not_genotyped)

    make_mock_wgs(outprefix, wgs_prefix, wgs_snps, origin_map_files.release_2024_now, new_sample_ids; 
        n_wgs_individuals=n_wgs_individuals,
    )
    common_snps_set = wgs_snps[1:n_common_snps]
    distinct_snps_sets = (;(mapfile_id => rand(setdiff(mapfile[:, 2], snp_intersection), n_distinct_snps)
        for (mapfile_id, mapfile) in zip(keys(origin_map_files), origin_map_files))...)
    
    snps_indices = Dict()
    for (mapfile_id, mapfile) in zip(keys(origin_map_files), origin_map_files)
        retained_snps = vcat(common_snps_set, distinct_snps_sets[mapfile_id])
        mapfile_snps_indices = sort(indexin(retained_snps, mapfile[:, 2]))
        subset_mapfile = mapfile[mapfile_snps_indices, :]
        snps_indices[mapfile_id] = mapfile_snps_indices
        write_map(outfile_prefix(outprefix, mapfile_id), subset_mapfile)
    end
    return snps_indices
end

function get_sample_ids(ped_lines) 
    map(ped_lines) do line
        sample_id_range = findfirst(r"\t.*?\t", line)
        line[sample_id_range.start+1:sample_id_range.stop-1]
    end
end

function get_genotype_description_start(line)
    tab_count = 0
    idx = 1
    while true
        try
            if line[idx] == '\t'
                tab_count += 1
            end
            if tab_count == 6
                return idx + 1
            end
        catch
            @warn("Skipping ill formed character, this should be solved by Dominique in the future.")
        end
        
        idx += 1
    end
end

function set_alleles!(alleles_A, alleles_B, ped_lines, variants_indices)
    for (col_id, line) in enumerate(ped_lines)
        first_variant_pos = get_genotype_description_start(line)
        for (row_id, variant_id) in enumerate(variants_indices)
            variant_pos = first_variant_pos + (variant_id - 1) * 4
            alleles_A[row_id, col_id] = line[variant_pos]
            alleles_B[row_id, col_id] = line[variant_pos+2]
        end
    end
end

function get_alleles(ped_lines, variants_indices)
    n_samples = size(ped_lines, 1)
    n_variants = size(variants_indices, 1)
    alleles_A = Matrix{Char}(undef, n_variants, n_samples)
    alleles_B = Matrix{Char}(undef, n_variants, n_samples)
    set_alleles!(alleles_A, alleles_B, ped_lines, variants_indices)
    return permutedims(alleles_A), permutedims(alleles_B)
end

function resample_variants(alleles_A, alleles_B)
    n_samples, n_variants = size(alleles_A)
    resampled_alleles_A = similar(alleles_A)
    resampled_alleles_B = similar(alleles_B)
    # Each variant is shuffled independently
    for variant_id in 1:n_variants
        perm = randperm(n_samples)
        resampled_alleles_A[:, variant_id] = alleles_A[perm, variant_id]
        resampled_alleles_B[:, variant_id] = alleles_B[perm, variant_id]
    end
    return permutedims(resampled_alleles_A), permutedims(resampled_alleles_B)
end

function make_new_ped_lines(
    ped_lines, resampled_alleles_A, 
    resampled_alleles_B, 
    new_sample_ids; 
    n_arrays_individuals=1000
    )
    new_ped_lines = String[]
    idx = 1
    while length(new_ped_lines) < n_arrays_individuals
        line = ped_lines[idx]
        # Get the line up to the genotypes
        line_prefix = line[1:get_genotype_description_start(line)-2]
        prefix_elements = split(line_prefix, '\t')
        # Substitute the new sample id
        prefix_elements[2] = new_sample_ids[idx]
        # Replace genotypes with resampled
        new_line = string(
            join(prefix_elements, '\t'),
            join('\t' .* resampled_alleles_A[:, idx] .* '\t' .* resampled_alleles_B[:, idx])
        )
        push!(new_ped_lines, new_line)
        idx += 1
    end
    return new_ped_lines
end

function make_mock_ped_files(
    release_r8, 
    release_2021_2023, 
    release_2024_now,
    new_sample_ids,
    variants_indices; 
    outprefix = "mock",
    n_arrays_individuals = 1000, 
    )
    origin_ped_files = (
        release_r8 = read_ped(release_r8),
        release_2021_2023 = read_ped(release_2021_2023),
        release_2024_now = read_ped(release_2024_now)
    );

    # Filter and Resample Variants
    for (file_id, ped_lines) in zip(keys(origin_ped_files), origin_ped_files)
        alleles_A, alleles_B = get_alleles(ped_lines, variants_indices[file_id])
        resampled_alleles_A, resampled_alleles_B = resample_variants(alleles_A, alleles_B)
        new_ped_lines = make_new_ped_lines(
            ped_lines, 
            resampled_alleles_A, 
            resampled_alleles_B, 
            new_sample_ids[file_id]; 
            n_arrays_individuals=n_arrays_individuals
        )
        write_ped(outfile_prefix(outprefix, file_id), new_ped_lines)
    end
end


function mock_genetic_data(
    release_r8,
    release_2021_2023, 
    release_2024_now,
    wgs_prefix,
    new_sample_ids;
    n_snps_wgs_not_genotyped=n_snps_wgs_not_genotyped,
    n_wgs_individuals=n_wgs_individuals,
    outprefix="mock",
    n_common_snps=100, 
    n_distinct_snps=10,
    n_arrays_individuals=1000,
    verbosity=1
    )

    # Make map files and wgs
    verbosity > 0 && @info("Creating mock map files.")
    variants_indices = make_mock_map_files_and_wgs(
        release_r8, 
        release_2021_2023, 
        release_2024_now,
        wgs_prefix,
        new_sample_ids; 
        outprefix=outprefix,
        n_common_snps=n_common_snps,
        n_distinct_snps=n_distinct_snps,
        n_wgs_individuals=n_wgs_individuals,
        n_snps_wgs_not_genotyped=n_snps_wgs_not_genotyped
    )
    # Make ped files
    verbosity > 0 && @info("Creating mock ped files.")
    make_mock_ped_files(
        release_r8, 
        release_2021_2023, 
        release_2024_now, 
        new_sample_ids, 
        variants_indices; 
        outprefix = outprefix,
        n_arrays_individuals = n_arrays_individuals
    )
end

function mock_covariates(covariates_path, new_sample_ids; outprefix="mock")
    covariates = CSV.read(covariates_path, DataFrame)
    id_pairs = collect(sample_id_map)
    sample_id_map = DataFrame([first.(id_pairs), last.(id_pairs)], ["OLD_ID", "NEW_ID"])
    covariates = innerjoin(covariates, sample_id_map, on=:genotype_file_id => :OLD_ID)
    select!(covariates,
        :NEW_ID => :genotype_file_id, 
        :age_years, 
        :sex, 
        :case_or_control, 
        :cohort, 
        :severe_cohort_primary_diagnosis, 
        :isaric_cohort_max_severity_score
    )
    CSV.write(string(outprefix, ".covariates.csv"), covariates)
end

function mock_data(release_r8,
    release_2021_2023, 
    release_2024_now,
    covariates_path,
    wgs_prefix;
    n_snps_wgs_not_genotyped = 30,
    n_wgs_individuals = 10,
    outprefix="mock",
    n_common_snps=100, 
    n_distinct_snps=10,
    n_arrays_individuals=1000,
    rng=123,
    verbosity=1)
    # Fix random seed
    Random.seed!(rng)
    # New sample ids
    new_sample_ids_ = string.("odap", 1:3*n_arrays_individuals+n_wgs_individuals)
    new_sample_ids = (
        release_r8 = new_sample_ids_[1:n_arrays_individuals],
        release_2021_2023 = new_sample_ids_[n_arrays_individuals+1:2*n_arrays_individuals],
        release_2024_now = new_sample_ids_[2*n_arrays_individuals+1:3*n_arrays_individuals],
        wgs = new_sample_ids_[3*n_arrays_individuals+1:end]
    )
    # Mock genetic data
    verbosity > 0 && @info("Mocking genotyping arrays.")
    mock_genetic_data(
        release_r8,
        release_2021_2023, 
        release_2024_now,
        wgs_prefix,
        new_sample_ids;
        n_snps_wgs_not_genotyped=n_snps_wgs_not_genotyped,
        n_wgs_individuals=n_wgs_individuals,
        outprefix=outprefix,
        n_common_snps=n_common_snps, 
        n_distinct_snps=n_distinct_snps,
        n_arrays_individuals=n_arrays_individuals,
        verbosity=verbosity-1
        )
    # Mock covariates
    verbosity > 0 && @info("Mocking covariates.")
    # mock_covariates(covariates_path, new_sample_ids; outprefix=outprefix)
end