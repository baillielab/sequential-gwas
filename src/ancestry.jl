function estimate_ancestry(genotypes_prefix, pedigree_file; output="ancestry.csv")
    # Write known ancestries to file
    fam = SequentialGWAS.read_fam(string(genotypes_prefix, ".fam"))
    pedigree = CSV.read(pedigree_file, DataFrame, select=[:SampleID, :Superpopulation])
    leftjoin!(fam, pedigree, on=:IID => :SampleID)
    replace!(fam.Superpopulation, missing => "-")
    open(string(genotypes_prefix, ".pop"), "w") do io
        for pop in fam.Superpopulation
            println(io, pop)
        end
    end
    # Run admixture
    K = length(unique(fam.Superpopulation)) - 1
    J = max(1, nthreads() - 1)
    run(`admixture $(genotypes_prefix).bed $K --supervised -j$J -s 123`)
    # Read ancestry indices (admixture does not make it clear which ancestry is which Q column)
    Q = readdlm(string(genotypes_prefix, ".$K.Q"))
    fam.SuperpopulationId = map(x -> x.I[2], argmax(Q, dims=2)[:, 1])
    # Map ancestry indices to populations using 1000 GP
    kgp = filter(row -> row.Superpopulation != "-", fam)
    id_to_pop = Dict()
    for (pop, id) in zip(kgp.Superpopulation, kgp.SuperpopulationId)
        if !haskey(id_to_pop, id)
            id_to_pop[id] = pop
        else
            @assert id_to_pop[id] == pop "Inconsistent assignment of populations"
        end
    end
    # Assign ancestries to all samples in our dataset and write to disk
    ancestry = filter(row -> row.Superpopulation == "-", fam)
    ancestry.Superpopulation = [id_to_pop[id] for id in ancestry.SuperpopulationId]
    CSV.write(
        output, 
        select(ancestry, :FID, :IID, :Superpopulation), 
    )
end