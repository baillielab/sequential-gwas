function tag_variant_id_missing_from_gwas!(pvar, gwas_results)
    # Map GWAS variants to their alleles
    gwas_ids_to_alleles = Dict(row.ID => Set([row.ALLELE0, row.ALLELE1]) for row in Tables.namedtupleiterator(select(gwas_results, [:ID, :ALLELE0, :ALLELE1])))
    # Update variants IDs in PVAR if they are not in the GWAS variants
    pvar.ID = map(Tables.namedtupleiterator(pvar)) do row
        if haskey(gwas_ids_to_alleles, row.ID)
            gwas_alleles = gwas_ids_to_alleles[row.ID]
            pvar_alleles = Set([row.ALT, row.REF])
            if gwas_alleles != pvar_alleles
                string(row.ID, ".unmatched_alleles")
            else
                row.ID
            end
        else
            string(row.ID, ".not_in_gwas")
        end
    end
end

function write_new_pgen_from_gwas_results(input_pgen_prefix, output_pgen_prefix, pvar, sample_file)
    # Write PVAr with tagged variants
    tmpdir = mktempdir()
    tagged_pvar_file = joinpath(tmpdir, "updated_variants.pvar")
    CSV.write(tagged_pvar_file, pvar; delim='\t')
    # Write variants to exclude to file
    exclude_list = filter(x -> endswith(x, ".not_in_gwas") || endswith(x, ".unmatched_alleles"), pvar.ID)
    exclude_file = joinpath(tmpdir, "exclude_variants.txt")
    open(exclude_file, "w") do io
        for id in exclude_list
            println(io, id)
        end
    end
    # Run PLINK filtering cmd
    run(`plink2 \
        --pgen $input_pgen_prefix.pgen \
        --psam $input_pgen_prefix.psam \
        --pvar $tagged_pvar_file \
        --make-pgen --out $output_pgen_prefix \
        --exclude $exclude_file \
        --keep $sample_file
    `)
end

function write_significant_clumps(pgen_prefix, gwas_results_file;
    min_sig_clump_size = 3,
    output = "clumps.sig.tsv",
    lead_pvalue = 5e-8,
    p2_pvalue = 5e-5,
    r2_threshold = 0.5,
    clump_kb = 250,
    clump_id_field = "ID",
    clump_pval_field = "LOG10P",
    allele_1_field = "ALLELE_1"
    )
    tmpdir = mktempdir()
    output_clump_prefix = joinpath(tmpdir, "clumps")
    run(`plink2 --pfile $pgen_prefix \
        --clump $gwas_results_file \
        --clump-p1 $lead_pvalue \
        --clump-p2 $p2_pvalue \
        --clump-r2 $r2_threshold \
        --clump-kb $clump_kb \
        --clump-id-field $clump_id_field \
        --clump-log10 \
        --clump-p-field $clump_pval_field \
        --clump-a1-field $allele_1_field \
        --out $output_clump_prefix
    `)

    clumps = CSV.read(output_clump_prefix * ".clumps", DataFrame; delim="\t")
    sig_clumps = filter(
        :SP2 => x -> x !== "." && length(split(x, ",")) >= min_sig_clump_size, 
        clumps
    )
    CSV.write(output, sig_clumps)
    return sig_clumps
end

function genotypes_from_pgen(pgen_prefix, ld_variants, sample_list)
    pgen = Pgen(string(pgen_prefix, ".pgen"))

    sample_idx = indexin(sample_list, pgen.psam_df.IID)
    nsamples = length(sample_idx)

    idx_inf = findfirst(==(ld_variants.ID_B[1]), pgen.pvar_df.ID)
    idx_sup = findfirst(==(ld_variants.ID_B[end]), pgen.pvar_df.ID)

    X = Matrix{Float64}(undef, nsamples, idx_sup - idx_inf + 1)
    g = Vector{UInt8}(undef, nsamples)
    g_ld = similar(g)
    col_id = 1
    for variant in iterator(pgen)
        if idx_inf <= variant.index <= idx_sup
            get_genotypes!(g, pgen, variant; ldbuf=g_ld)
            v_rt = variant.record_type & 0x07
            if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
                g_ld .= g
            end
            variant_genotypes = g[sample_idx]
            μ = mean(filter(!=(0x03), variant_genotypes))
            X[:, col_id] .= [g_ === 0x03 ? μ : g_ for g_ in variant_genotypes]
            col_id += 1
        end
    end
    return X, pgen.pvar_df[idx_inf:idx_sup, :]
end

function dosages_from_pgen()
    d = Vector{Float32}(undef, n_samples(p))
    g = Vector{UInt8}(undef, n_samples(p))
    g_ld = similar(g)
    for v in v_iter
        alt_allele_dosage!(d, g, p, v; genoldbuf=g_ld)
        v_rt = v.record_type & 0x07
        if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
            g_ld .= g
        end
        
        # do someting with dosage values in `d`...
    end
end

function get_phenotype(covariates_file, sample_list, phenotype)
    covariates = CSV.read(covariates_file, DataFrame; delim='\t', select=["IID", phenotype])
    sample_idx = indexin(sample_list, covariates.IID)
    return covariates[sample_idx, phenotype]
end

function susie_finemap(X, y; n_causal=10)
    @rput X
    @rput y
    @rput n_causal
    R"""
    library(susieR)
    fitted = susie(X, y, L=n_causal)
    """
    return @rget fitted
end

function get_ld_variants(variant_id, pgen_prefix; ld_window_kb=1000, ld_window_r2=0.1)
    tmpdir = mktempdir()
    out_prefix = joinpath(tmpdir, "$variant_id.LD")
    run(`plink2 \
        --pfile $pgen_prefix \
        --r2-phased \
        --ld-snp $variant_id \
        --ld-window-kb $ld_window_kb \
        --ld-window 99999 \
        --ld-window-r2 $ld_window_r2 \
        --out $out_prefix
    `)
    return CSV.read(out_prefix * ".vcor", DataFrame, delim='\t')
end

update_credible_sets!(cs_vector, variant_idx::Int, cs::Symbol) = cs_vector[variant_idx] = parse(Int, replace(string(cs), "L" => ""))

function update_credible_sets!(cs_vector, variant_idxs::AbstractVector, cs::Symbol)
    for variant_idx in variant_idxs
        update_credible_sets!(cs_vector, variant_idx, cs)
    end
end

function update_credible_sets!(cs_vector, susie_results)
    for (cs, variant_idxs) in susie_results[:sets][:cs]
        update_credible_sets!(cs_vector, variant_idxs, cs)
    end
end

"""
    finemap_significant_regions(
        gwas_results_file,
        pgen_prefix,
        covariates_file,
        sample_file;
        output_prefix = "clumps.sig.tsv",
        min_sig_clump_size = 3,
        lead_pvalue = 5e-8,
        p2_pvalue = 5e-5,
        r2_threshold = 0.5,
        clump_kb = 250,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1",
        ld_window_kb=1000,
        ld_window_r2=0.1,
        n_causal = 10
        )

This function performs fine-mapping of significant regions identified from GWAS results.

1. Identifies variants from the PGEN fileset that have gone through GWAS (some may not due to MAF/MAC/...)
2. Identify clumps of significant variants from the GWAS results to tag independent association regions
3. Finemap these regions using SuSiE

# Arguments
- `gwas_results_file::String`: Path to the GWAS results file (TSV format).
- `pgen_prefix::String`: Prefix for the PGEN fileset (without .pgen extension).
- `covariates_file::String`: Path to the covariates file (TSV format).
- `sample_file::String`: Path to the sample IDs file used to generate the GWAS results.
- `output_prefix::String`: Prefix to output the significant clumps (TSV format).
- `min_sig_clump_size::Int`: Minimum number of variants in a clump to be considered significant.
- `lead_pvalue::Float64`: P-value threshold for lead variants in clump.
- `p2_pvalue::Float64`: Secondary p-value threshold for variants in clump.
- `r2_threshold::Float64`: LD r² threshold for clumping.
- `clump_kb::Int`: Distance in kb for clumping.
- `clump_id_field::String`: Column name in GWAS results for variant IDs.
- `clump_pval_field::String`: Column name in GWAS results for p-values.
- `allele_1_field::String`: Column name in GWAS results for allele 1.
- `ld_window_kb::Int`: Window size in kb for LD calculation around lead variant to define the fine mapping region.
- `ld_window_r2::Float64`: r² threshold for including variants in LD window to define the fine mapping region.
- `n_causal::Int`: Number of causal variants to assume in SuSiE fine-mapping.
"""
function finemap_significant_regions(
    gwas_results_file,
    pgen_prefix,
    covariates_file,
    sample_file;
    output_prefix = "finemapping_results",
    min_sig_clump_size = 3,
    lead_pvalue = 5e-8,
    p2_pvalue = 5e-5,
    r2_threshold = 0.5,
    clump_kb = 250,
    clump_id_field = "ID",
    clump_pval_field = "LOG10P",
    allele_1_field = "ALLELE_1",
    ld_window_kb=1000,
    ld_window_r2=0.1,
    n_causal = 10
    )
    _, _, group, phenotype, _ = split(gwas_results_file, ".")
    # Read GWAS results and PVAR files
    @info "Reading GWAS results and PVAR files"
    gwas_results = CSV.read(gwas_results_file, DataFrame)
    pvar = CSV.read(pgen_prefix * ".pvar", DataFrame; delim='\t', comment="##")
    # Tag variant IDs that are not in the GWAS results
    @info "Tagging variants not in GWAS results"
    tag_variant_id_missing_from_gwas!(pvar, gwas_results)

    # Make new PGEN fileset with variants from GWAS only
    @info "Writing new PGEN fileset with variants matched to GWAS results"
    tmpdir = mktempdir()
    gwas_matched_pgen_prefix = joinpath(tmpdir, "gwas_matched")
    write_new_pgen_from_gwas_results(pgen_prefix, gwas_matched_pgen_prefix, pvar, sample_file)

    # Find clumps
    @info "Finding clumps in GWAS-matched PGEN fileset"
    sig_clumps = write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = min_sig_clump_size,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = lead_pvalue,
        p2_pvalue = p2_pvalue,
        r2_threshold = r2_threshold,
        clump_kb = clump_kb,
        clump_id_field = clump_id_field,
        clump_pval_field = clump_pval_field,
        allele_1_field = allele_1_field
    )

    sample_list = getindex.(split.(readlines(sample_file), "\t"), 2)
    
    # Finemap each clump
    finemapping_results = []
    for clump in eachrow(sig_clumps)
        @info "Fine Mapping clump led by : $(clump.ID)"
        ld_variants = get_ld_variants(clump.ID, gwas_matched_pgen_prefix; ld_window_kb=ld_window_kb, ld_window_r2=ld_window_r2)
        X, variants_info = genotypes_from_pgen(gwas_matched_pgen_prefix, ld_variants, sample_list)
        y = get_phenotype(covariates_file, sample_list, phenotype)
        susie_results = susie_finemap(X, y; n_causal=n_causal)
        variants_info.PIP = susie_results[:pip]
        variants_info.CS = zeros(Int, size(X, 2))
        variants_info.CLUMP_ID = fill(clump.ID, size(X, 2))
        update_credible_sets!(variants_info.CS, susie_results)
        variants_info = innerjoin(
            variants_info, 
            gwas_results[!, [:ID, :ALLELE0, :ALLELE1, :A1FREQ, :N, :TEST, :BETA, :SE, :CHISQ, :LOG10P]], 
            on=:ID
        )
        push!(finemapping_results, variants_info)
    end
    CSV.write(string(output_prefix, ".tsv"), vcat(finemapping_results...))

    return 0
end

