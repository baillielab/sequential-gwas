function make_sample_batches(prefix)
    fam = SequentialGWAS.read_fam(string(prefix, ".fam"))
    return map(Iterators.partition(1:nrow(fam), samples_per_file)) do indices
        filename = string(prefix, ".samples_", indices[1], "_", indices[end], ".keep")
        CSV.write(filename, fam[indices, [:FID, :IID]], delim=' ', header=false)
        filename
    end
end

function get_chromosomes(prefix)
    bim = SequentialGWAS.read_bim(string(prefix, ".bim"))
    return unique(bim.CHR_CODE)
end

function split_bed_file_to_vcf(prefix)
    sample_batches = make_sample_batches(prefix)
    chromosomes = get_chromosomes(prefix)
    Threads.@threads for chr in chromosomes
        for samples_batch in sample_batches
            samples_ext = splitext(first(splitext(basename(samples_batch))))[end]
            temp_prefix = string(prefix, ".", chr, samples_ext)
            run(pipeline(
                Cmd([
                    "plink2",
                    "--bfile", prefix,
                    "--threads", "1",
                    "--output-chr", "chr26",
                    "--keep", samples_batch,
                    "--chr", string(chr),
                    "--recode", "vcf",
                    "--out", temp_prefix
                ]), 
                Cmd(
                    ["bgzip", 
                    string(temp_prefix, ".vcf")
                ]),
            ))
        end
    end
end

function impute(prefix, token_file; 
    password="abcde", 
    r2_filter=0.8, 
    submission_batchsize=3,
    refresh_rate=60, 
    samples_per_file=10_000
    )
    prefix = "/Users/olabayle/Dev/sequential-gwas/test/assets/gwas/genotypes/genotypes.arrays_wgs.aggregated"
    split_bed_file_to_vcf(prefix)
    token_file = "assets/topmed-api-token"
    token = read(token_file, String)
    dir, _prefix = splitdir(prefix)
    vcf_files = filter(x -> endswith(x, ".vcf.gz"), readdir(dir, join=true))
    (batchid, batch) = first(enumerate(Iterators.partition(vcf_files, batchsize)))
    for (batchid, batch) in enumerate(Iterators.partition(vcf_files, batchsize))
        # Submit
        cmd = Cmd([
                "curl", 
                "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver",
                "-X", "POST",
                "-H", "X-Auth-Token: $token",
                "-F", "job-name=$batchid",
                [string("-F files=@$(realpath(f))") for f in batch]...,
                "-F", "refpanel=apps@topmed-r3",
                "-F", "build=hg38",
                "-F", "phasing=eagle",
                "-F", "password=$password",
                "-F", "population=all",
                "-F", "r2Filter=$r2_filter",
                "-F", "meta=yes"
        ])
        job_details = JSON.parse(read(cmd, String))
        job_details["success"] == true || throw(error(job_details["message"]))
        job_id = job_details["id"]
        # Wait for response
        while true
            cmd = Cmd([
                "curl", 
                "-H", "X-Auth-Token: $token", 
                string("https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/$job_id/status")
            ])
            status = JSON.parse(read(cmd, String))
            state = status["state"]
            if state == 5 || state == 6
                throw(error("Job $batchid failed."))
            elseif state == 4
                break
            else
                sleep(rate)
            end
        end
        # Download results
    end
end