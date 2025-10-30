outfile_prefix(prefix, file_id) = string(prefix, ".", file_id)

function make_variants_lists(origin_map_files;
    n_common_snps=100,
    n_distinct_snps=10,
    n_snps_wgs_not_genotyped=30
    )
    wgs_variants = []
    genotyping_arrays_common_variants = []
    genotyping_arrays_unique_variants = Dict(
        :release_r8 => [],
        :release_2021_2023 => [],
        :release_2024_now => []
    )
    for chr in 1:22
        release_r8_variants = filter(:CHR => ==(chr), origin_map_files.release_r8).ID
        release_2021_2023_variants = filter(:CHR => ==(chr), origin_map_files.release_2021_2023).ID
        release_2024_now_variants = filter(:CHR => ==(chr), origin_map_files.release_2024_now).ID
        common_variants = intersect(release_r8_variants, release_2021_2023_variants, release_2024_now_variants)
        # Add chr variants to wgs_variants
        chr_wgs_variants = rand(common_variants, n_common_snps+n_snps_wgs_not_genotyped)
        append!(wgs_variants, chr_wgs_variants)
        # Add chr variants to genotyping_arrays_common_variants
        append!(genotyping_arrays_common_variants, chr_wgs_variants[1:n_common_snps])
        # Add unique variants per release
        genotyping_arrays_unique_variants[:release_r8] = rand(setdiff(release_r8_variants, common_variants), n_distinct_snps)
        genotyping_arrays_unique_variants[:release_2021_2023] = rand(setdiff(release_2021_2023_variants, common_variants), n_distinct_snps)
        genotyping_arrays_unique_variants[:release_2024_now] = rand(setdiff(release_2024_now_variants, common_variants), n_distinct_snps)
    end
    return wgs_variants, genotyping_arrays_common_variants, genotyping_arrays_unique_variants
end

function make_mock_wgs(outprefix, wgs_prefix, wgs_snps, release_2024_now_map, new_sample_ids; 
    n_wgs_individuals=10,
    )
    # Make regions file to extract
    map_df = DataFrame(release_2024_now_map, ["CHR", "ID", "POS", "BP"])
    wgs_snps_df = DataFrame(ID=wgs_snps)
    merged = innerjoin(DataFrames.select(map_df, [:CHR, :ID, :BP]), wgs_snps_df, on=:ID)
    merged.CHR .= string.("chr", merged.CHR)
    merged.BP_MINUS_ONE = merged.BP .- 1
    merged.BP_PLUS_ONE = merged.BP .+ 1
    tmpdir = mktempdir()
    regions_file = joinpath(tmpdir, "variants.tsv")
    CSV.write(regions_file, DataFrames.select(merged, [:CHR, :BP_MINUS_ONE, :BP_PLUS_ONE]), delim='\t', header=false)
    # List gvcf files
    gvcf_dir = dirname(wgs_prefix)
    gvcf_files = filter(endswith("gvcf.gz"), readdir(gvcf_dir))
    gvcf_files = rand(gvcf_files, n_wgs_individuals)
    # Extract regions from gvcf files
    input_gvcf = gvcf_files[2]
    new_sample_id = new_sample_ids.wgs[2]
    for (input_gvcf, new_sample_id) in zip(gvcf_files, new_sample_ids.wgs)
        input_gvcf_filepath = joinpath(gvcf_dir, input_gvcf)
        output_gvcf_filepath_temp = joinpath(tmpdir, string("temp_gvcf_", new_sample_id, ".gvcf.gz"))
        output_gvcf_filepath = string(outprefix, ".", new_sample_id, ".gvcf.gz")
        new_sample_file = joinpath(tmpdir, string("new_sample_", new_sample_id ,".txt"))
        write(new_sample_file, new_sample_id)
        run(`bcftools view -O z -R $regions_file -o $output_gvcf_filepath_temp $input_gvcf_filepath`)
        run(`bcftools reheader -s $new_sample_file -o $output_gvcf_filepath $output_gvcf_filepath_temp`)
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
        release_r8 = read_map(string(release_r8, ".map")),
        release_2021_2023 = read_map(string(release_2021_2023, ".map")),
        release_2024_now = read_map(string(release_2024_now, ".map"))
    )

    wgs_variants, genotyping_arrays_common_variants, genotyping_arrays_unique_variants = make_variants_lists(
        origin_map_files; 
        n_common_snps=n_common_snps,
        n_distinct_snps=n_distinct_snps,
        n_snps_wgs_not_genotyped=n_snps_wgs_not_genotyped
    )

    make_mock_wgs(outprefix, wgs_prefix, wgs_variants, origin_map_files.release_2024_now, new_sample_ids; 
        n_wgs_individuals=n_wgs_individuals,
    )

    snps_indices = Dict()
    for (mapfile_id, mapfile) in zip(keys(origin_map_files), origin_map_files)
        retained_snps = vcat(genotyping_arrays_common_variants, genotyping_arrays_unique_variants[mapfile_id])
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
    verbosity=1
    )
    origin_ped_files = (
        release_r8 = read_ped(release_r8),
        release_2021_2023 = read_ped(release_2021_2023),
        release_2024_now = read_ped(release_2024_now)
    );

    # Filter and Resample Variants
    for (file_id, ped_lines) in zip(keys(origin_ped_files), origin_ped_files)
        verbosity > 0 && @info(string("creating mock ped file: ", file_id))
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
        n_arrays_individuals = n_arrays_individuals,
        verbosity=verbosity
    )
end

function mock_covariates(covariates_path, new_sample_ids; outprefix="mock")
    covariates = CSV.read(covariates_path, DataFrame)
    new_ids = vcat(new_sample_ids...)
    covariates_sample = shuffle(covariates)[1:length(new_ids), :]
    covariates_sample.genotype_file_id = new_ids
    select!(covariates_sample,
        :genotype_file_id, 
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
    mock_covariates(covariates_path, new_sample_ids; outprefix=outprefix)
end