
function write_new_bim_bed_fam(input_prefix, shared_variants_file, output_prefix)
    run(Cmd([
        "plink2",
        "--threads", string(nthreads()),
        "--bfile", input_prefix,
        "--max-alleles", "2",
        "--output-chr", "chr26",
        "--extract", shared_variants_file,
        "--make-bed",
        "--out", output_prefix
    ]))
end

function vcf_to_bed(vcf_file, output_prefix)
        # Converts to plink format
    run(Cmd(["plink2",
            "--threads", string(nthreads()),
            "--vcf", vcf_file,
            "--output-chr", "chr26",
            "--max-alleles", "2",
            "--set-all-var-ids", "@:#",
            "--make-bed",
            "--out", output_prefix
    ]))

end

function gatk_genotyping(vcf_file, shared_variants_gatk, reference_genome, output_vcf)
    run(Cmd(["gatk", "GenotypeGVCFs",
    "-R", reference_genome,
    "-V", vcf_file,
    "--intervals", shared_variants_gatk,
    "--interval-padding", "10",
    "--force-output-intervals", shared_variants_gatk,
    "-O", output_vcf
    ]))
end

function allele_map_from_shared_variants(shared_variants_file)
    variants_map = Dict()
    for line in readlines(shared_variants_file)
        chr, pos, ref, alt = split(line, ":")
        variants_map[string(chr, ":", pos)] = (
            variant_id = line,
            allele_map = Dict(ref => alt, alt => ref)
        )
    end
    return variants_map
end

function update_bim_with_mapped_alleles(tmp_prefix, shared_variants_file)
    bim_file = string(tmp_prefix, ".bim")
    # Convert VCF to a temporary PLINK bed
    variants_map = allele_map_from_shared_variants(shared_variants_file)
    # Fill the missing allele and update the variant id in the bim file
    bim = SequentialGWAS.read_bim(bim_file)
    bim.VARIANT_ID = convert(Vector{String}, bim.VARIANT_ID)
    bim.ALLELE_1 = convert(Vector{String}, bim.ALLELE_1)
    for row in eachrow(bim)
        ## Check if the variant is in the shared variants (we also genotype around the variants so it could not be)
        if haskey(variants_map, row.VARIANT_ID)
            info = variants_map[row.VARIANT_ID]
            if row.ALLELE_1 == "."
                try
                    row.ALLELE_1 = info.allele_map[row.ALLELE_2]
                catch
                    ## If the individual has a non ref/alt variant from the genotyping arrays (i.e. rare) we will not integrate it
                    return row.VARIANT_ID
                end
            end
            row.VARIANT_ID = info.variant_id
        end
    end
    # Write the variants that are not in the shared variants
    CSV.write(bim_file, bim, header=false, delim="\t")
    return "OK"
end

"""

    genotype_gvcf(gvcf_file, 
        shared_variants_plink, 
        shared_variants_gatk, 
        reference_genome; 
        output_prefix="output"
    )

Genotype a GVCF file using GATK and PLINK as follow. 

- First attempts to genotype the GVCF file using GATK (this can fail for some reason in which case we just skip the individual).
- Then converts the VCF to a PLINK bed format.
- Updates the bim file with the mapped alleles from the shared variants (thi can also fail if some variants are not any of the known ref/alt in which case we also skip the individual).
- Writes the new bim, bed and fam files.
"""
function genotype_gvcf(gvcf_file, 
    shared_variants_plink, 
    shared_variants_gatk, 
    reference_genome; 
    output_prefix="output"
    )
    # All temporary files will be stored in a temporary directory
    tmp_dir = mkdir("temp")
    tmp_prefix = joinpath(tmp_dir, "tmp")
    tmp_vcf = string(tmp_prefix, ".vcf.gz")
    # Genotype GVCF using GATK
    gatk_genotyping(gvcf_file, shared_variants_gatk, reference_genome, tmp_vcf)
    # Convert VCF to a temporary PLINK bed
    vcf_to_bed(tmp_vcf, tmp_prefix)
    # Update bim file with the mapped alleles
    status = update_bim_with_mapped_alleles(tmp_prefix, shared_variants_plink)
    # If the individual has a non ref/alt variant from the genotyping arrays (i.e. rare) we will not integrate it
    if status != "OK"
        write(string(output_prefix, ".txt"), status)
    else
        write_new_bim_bed_fam(tmp_prefix, shared_variants_plink, output_prefix)
    end
    # Clean
    rm(tmp_dir; recursive=true)
    return 0
end