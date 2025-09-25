module PostGWAS

using CSV
using DataFrames
using Parquet
using TMLE
using Combinatorics
using RCall

function extract_variant_info(variant_id)
    chr, pos, ref, alt, assembly = split(variant_id, "_")
    return (CHR=replace(chr, "chr" => ""), POS=parse(Int, pos), REF=ref, ALT=alt)
end

join_variant_info(chr_col, pos_col, ref_col, alt_col) = 
    string.(string.(chr_col), ":", string.(pos_col), ":", string.(ref_col), ":", string.(alt_col))

function add_universal_ID!(df; cols=[:CHR, :POS, :REF, :ALT])
    transform!(df, 
        cols => join_variant_info => :UNIV_ID)
    return df
end

function get_harmonized_eqtl_data!(parquet_file)
    eqtl_data = read_parquet(parquet_file) |> DataFrame
    tissue = split(basename(parquet_file), ".")[1]
    transform!(eqtl_data, :variant_id => ByRow(extract_variant_info) => AsTable)
    add_universal_ID!(eqtl_data)
    eqtl_data.TISSUE .= tissue
    return eqtl_data
end

function get_harmonized_finemapping_results!(finemapping_file)
    finemapping_results = CSV.read(finemapping_file, DataFrame, delim="\t")
    rename!(finemapping_results, "#CHROM" => "CHR")
    add_universal_ID!(finemapping_results)
    return finemapping_results
end

function map_variants_to_GTEx_genes(variants_df, parquet_files)
    variant_target_pairs = []
    # parquet_file = "../assets/GTEx_Analysis_v10_eQTL_updated/Artery_Tibial.v10.eQTLs.signif_pairs.parquet"
    for parquet_file in parquet_files
        sig_cis_eqtls = get_harmonized_eqtl_data!(parquet_file)
        merged_results = innerjoin(
            select(variants_df, :UNIV_ID),
            select(sig_cis_eqtls, 
                :UNIV_ID,
                :TISSUE,
                :gene_id => :GENE_ID, 
                :pval_nominal => :EQTL_PVAL, 
                :slope => :EQTL_BETA, 
                :slope_se => :EQTL_BETA_SE
            ), 
            on =:UNIV_ID
        )
        push!(variant_target_pairs, merged_results)
    end
    return vcat(variant_target_pairs...)
end

function download_GTex_QTLs(;file="GTEx_Analysis_v10_eQTL.tar", output_dir="assets/")
    tar_file = joinpath(output_dir, file)
    output_file = joinpath(output_dir, replace(file, ".tar" => "_updated"))
    if !isfile(output_file)
        Downloads.download("https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/$file", tar_file)
        run(`tar -xf $(tar_file) -C $(output_dir)`)
        rm(tar_file)
    end
end

function download_all_GTEx_QTLs_v10(;output_dir="assets/")
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    download_GTex_QTLs(file="GTEx_Analysis_v10_eQTL.tar", output_dir=output_dir)
    download_GTex_QTLs(file="GTEx_Analysis_v10_sQTL.tar", output_dir=output_dir)
end

function julia_main(gwas_file, finemapping_file, gtex_dir; output_prefix="post_gwas")
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "post_gwas")
    phenotype = "SEVERE_COVID_19"
    gwas_file = "test/assets/results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"
    finemapping_file = "test/assets/results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"
    gtex_dir = "../assets/GTEx_Analysis_v10_eQTL_updated"
    confounders = ["PC1", "PC2", "PC3", "PC4", "PC5"]
    outcome_extra_covariates = ["AGE", "SEX"]
    positivity_constraint = 0.01

    parquet_files = filter(endswith(".parquet"), readdir(gtex_dir, join=true))
    gwas_results = CSV.read(gwas_file, DataFrame, delim="\t")
    finemapping_results = get_harmonized_finemapping_results!(finemapping_file)
    filter!(x -> x.CS !== missing, finemapping_results)

    # Map finemapped variants to the genes they are eQTLs for
    variant_target_pairs = map_variants_to_GTEx_genes(finemapping_results, parquet_files)
    CSV.write(string(output_prefix, ".finemapped_variant_egene_pairs.tsv"), variant_target_pairs; delim="\t", header=true)

    # Get interaction candidates
    egenes_with_more_than_one_variant = filter(
        :UNIV_IDs => x -> length(split(x, ",")) > 1,
        combine(
            groupby(variant_target_pairs, [:GENE_ID, :TISSUE]), 
            :UNIV_ID => (x -> join(x, ",")) => :UNIV_IDs
        )
    )
    egenes_with_more_than_one_variant.ESTIMANDS = map(egenes_with_more_than_one_variant.UNIV_IDs) do variant_str
        variants = split(variant_str, ",")
        map(combinations(variants, 2)) do variants_comb
            factorialEstimand(
                AIE,
                variants_comb,
                phenotype,
                dataset=dataset,
                confounders=confounders,
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint
            )
        end
    end
end

end # module PostGWAS
