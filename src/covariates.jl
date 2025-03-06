function combine_covariates(
    covariates_file, 
    ancestry_file, 
    pcs_file; 
    output="covariates.merged.csv"
    )
    covariates = CSV.read(covariates_file, DataFrame)
    rename!(covariates, :genotype_file_id => :IID)
    ancestries = CSV.read(ancestry_file, DataFrame, select=[:IID, :Superpopulation])
    pcs = CSV.read(pcs_file, DataFrame, drop=["#FID"])
    merged = innerjoin(
        innerjoin(
            covariates, 
            ancestries, 
            on=:IID
        ),
        pcs,
        on=:IID
    )
    CSV.write(output, merged)
    return merged
end