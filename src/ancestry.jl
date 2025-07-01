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

function estimate_ancestry(genotypes_prefix, pedigree_file; output="ancestry.csv", threshold=0.8)
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