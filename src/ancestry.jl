function add_pedigree_info_and_write_to_file!(fam, pedigree_file, genotypes_prefix)
    pedigree = CSV.read(pedigree_file, DataFrame, select=[:SampleID, :Superpopulation])
    pedigree = Dict(zip(pedigree.SampleID, pedigree.Superpopulation))
    open(string(genotypes_prefix, ".pop"), "w") do io
        fam.Superpopulation = map(fam.IID) do sample_id
            pop = get(pedigree, sample_id, "-")
            println(io, pop)
            pop
        end
    end
end

function ancestry_from_fam(fam, index_to_pop; threshold=0.8)
    ancestry = filter(row -> row.Superpopulation == "-", fam)
    ancestry.Superpopulation = map(zip(ancestry.MostLikelyAncestryIndex, ancestry.MostLikelyAncestryProba)) do (cart_index, proba)
        if proba > threshold
            index_to_pop[cart_index[2]]
        else
            "ADMIXED"
        end
    end
    return ancestry
end


function estimate_ancestry(genotypes_prefix, pedigree_file; mode="admixture", output="ancestry.csv", threshold=0.8)
    if mode == "admixture"
        return admixture_ancestry_estimation(genotypes_prefix, pedigree_file; output=output, threshold=threshold)
    elseif mode == "scope"
        return scope_ancestry_estimation(genotypes_prefix, pedigree_file; output=output, threshold=threshold)
    else
        error("Unknown mode: $mode. Use 'admixture' or 'scope'.")
    end
end

function admixture_ancestry_estimation(genotypes_prefix, pedigree_file; output="ancestry.csv", threshold=0.8)
    # Write known ancestries to file
    fam = GenomiccWorkflows.read_fam(string(genotypes_prefix, ".fam"))
    GenomiccWorkflows.add_pedigree_info_and_write_to_file!(fam, pedigree_file, genotypes_prefix)
    # Run admixture
    K = length(unique(fam.Superpopulation)) - 1
    J = max(1, nthreads() - 1)
    run(`admixture $(genotypes_prefix).bed $K --supervised -j$J -s 123`)
    # Read ancestry indices (admixture does not make it clear which ancestry is which Q column)
    Q = readdlm(string(basename(genotypes_prefix), ".$K.Q"))
    fam.MostLikelyAncestryIndex = argmax(Q, dims=2)[:, 1]
    fam.MostLikelyAncestryProba = Q[fam.MostLikelyAncestryIndex]
    # Map ancestry indices to populations using 1000 GP
    kgp = filter(row -> row.Superpopulation != "-", fam)
    index_to_pop = Dict()
    for (pop, cart_index) in zip(kgp.Superpopulation, kgp.MostLikelyAncestryIndex)
        col_index = cart_index[2]
        if !haskey(index_to_pop, col_index)
            index_to_pop[col_index] = pop
        else
            @assert index_to_pop[col_index] == pop "Inconsistent assignment of populations"
        end
    end
    for (colindex, Qcol) in enumerate(eachcol(Q))
        fam[!, index_to_pop[colindex]] = Qcol
    end
    # Assign ancestries to all samples in our dataset and write to disk
    ancestry = ancestry_from_fam(fam, index_to_pop; threshold=threshold)
    CSV.write(
        output, 
        DataFrames.select(ancestry, :FID, :IID, :Superpopulation, values(index_to_pop)...), 
    )
end


function read_scope_ancestry_estimates(n_indiv, scope_ancestry_file)
    Q_lines = strip.(readlines(scope_ancestry_file))
    Q = Matrix{Float64}(undef, n_indiv, 5)
    for (pop_index, line) in enumerate(Q_lines)
        striped_line = filter.(!=(""), split.(line, isspace))
        Q[:, pop_index] = parse.(Float64, striped_line)
    end
    return Q
end

function assign_scope_ancestry_estimates!(fam, Q; threshold=0.8, ordered_ancestries = ["AFR", "AMR", "EAS", "EUR", "SAS"])
    for (ancestry, estimate) in zip(ordered_ancestries, eachcol(Q))
        fam[!, ancestry] = estimate
    end
    # Assign most likely ancestry
    most_likely_ancestries_index = getindex.(argmax(Q, dims=2), 2)
    max_probas = maximum(Q, dims=2)
    most_likely_ancestries = ordered_ancestries[most_likely_ancestries_index]
    fam.Superpopulation = map(zip(max_probas[:, 1], most_likely_ancestries[:, 1])) do (max_proba, mostl_likely_ancestry)
        max_proba > threshold ? mostl_likely_ancestry : "ADMIXED"
    end
    return fam
end

function extract_population(input_prefix, individuals; tmpdir=mktempdir(), output_prefix="kgp")
    keep_file = joinpath(tmpdir, string(output_prefix, ".indiv.txt"))
    CSV.write(keep_file, individuals[!, [:FID, :IID]], header=false, delim="\t")
    full_output_prefix = joinpath(tmpdir, output_prefix)
    run(`plink --bfile $input_prefix --keep $keep_file --make-bed --out $full_output_prefix`)
    return full_output_prefix
end


function format_stratified_freqs(input_file, output_file)
    cluster_map = Dict("AFR" => "1", "AMR" => "2", "EAS" => "3", "EUR" => "4", "SAS" => "5")
    lines = strip.(readlines(input_file))
    open(output_file, "w") do io
        println(io, "CHR\tSNP\tCLST\tA1\tA2\tMAF\tMAC\tNCHROBS")
        for line in lines[2:end]
            striped_line = filter.(!=(""), split.(line, isspace))
            striped_line[3] = cluster_map[striped_line[3]]
            striped_line[7] = "N"
            striped_line[8] = "N"
            println(io, join(striped_line, "\t"))
        end
    end
    return output_file
end

function scope_ancestry_estimation(genotypes_prefix, kgp_pedigree_file; output="ancestry.csv", threshold=0.8)
    tmpdir = mktempdir()
    fam = GenomiccWorkflows.read_fam(string(genotypes_prefix, ".fam"))
    kgp_pedigrees = CSV.read(kgp_pedigree_file, DataFrame, select=[:SampleID, :Superpopulation])
    # Write KGP genotypes
    kgp_individuals = filter(x -> x.IID ∈ kgp_pedigrees.SampleID, fam)
    kgp_prefix = GenomiccWorkflows.extract_population(genotypes_prefix, kgp_individuals; tmpdir=tmpdir, output_prefix="kgp")
    # Write Other cohort genotypes
    other_individuals = filter(x -> x.IID ∉ kgp_individuals.IID, fam)
    other_prefix = GenomiccWorkflows.extract_population(genotypes_prefix, other_individuals; tmpdir=tmpdir, output_prefix="other")
    # Write KGP superpopulation 
    kgp_fam = GenomiccWorkflows.read_fam(string(kgp_prefix, ".fam"))
    pedigree = Dict(zip(kgp_pedigrees.SampleID, kgp_pedigrees.Superpopulation))
    kgp_fam.SuperPopulation = map(iid -> pedigree[iid], kgp_fam.IID)
    pop_file = joinpath(tmpdir, "kgp.pop.tsv")
    CSV.write(pop_file, kgp_fam[!, [:FID, :IID, :SuperPopulation]], header=false, delim="\t")
    # Compute stratified KGP frequencies with plink
    output_freqs_prefix = joinpath(tmpdir, "kgp")
    run(`plink --bfile $kgp_prefix --freq --within $pop_file --out $output_freqs_prefix`)
    # Adapt the poorly formated frequency file for SCOPE, see: https://github.com/sriramlab/SCOPE/blob/master/misc/simulations/generate_plink_frq.py
    scope_freq_file = format_stratified_freqs(
        string(output_freqs_prefix, ".frq.strat"), 
        joinpath(tmpdir, "kgp.SCOPE.frq")
    )
    # Run SCOPE
    output_prefix = joinpath(tmpdir, "scope.results")
    run(`scope -g $other_prefix -k 5 -seed 123 -freq $scope_freq_file -o $output_prefix`)
    # Read SCOPE results
    other_fam = GenomiccWorkflows.read_fam(string(other_prefix, ".fam"))
    Q = read_scope_ancestry_estimates(nrow(other_fam), string(output_prefix, "Qhat.txt"))
    # Assign ancestry estimates
    ordered_ancestries = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    assign_scope_ancestry_estimates!(other_fam, Q; threshold=threshold, ordered_ancestries=ordered_ancestries)
    # Write Output Table
    CSV.write(
        output, 
        DataFrames.select(other_fam, :FID, :IID, :Superpopulation, ordered_ancestries...), 
        header=true
    )
    # Clean up
    rm(tmpdir, recursive=true)
    
    return 0
end