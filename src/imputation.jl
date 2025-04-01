function make_sample_batches(prefix; output_dir=".", samples_per_file=5_000)
    fam = SequentialGWAS.read_fam(string(prefix, ".fam"))
    return map(Iterators.partition(1:nrow(fam), samples_per_file)) do indices
        filename = joinpath(output_dir, string("samples_", indices[1], "_", indices[end], ".keep"))
        CSV.write(filename, fam[indices, [:FID, :IID]], delim=' ', header=false)
        filename
    end
end

function get_chromosomes(prefix)
    bim = SequentialGWAS.read_bim(string(prefix, ".bim"))
    return unique(bim.CHR_CODE)
end

function get_vcf_files_channel(dir)
    vcf_files = filter(x -> endswith(x, ".vcf.gz"), readdir(dir, join=true))
    samples = [join(split(replace(basename(f), ".vcf.gz" => ""), "_")[3:4], "-") for f in vcf_files]
    chromosomes = [split(replace(basename(f), ".vcf.gz" => ""), "_")[1] for f in vcf_files]
    vcf_files_df = DataFrame(FILE=vcf_files, SAMPLES=samples, CHR=chromosomes)
    return Channel() do channel
        for (key, group) in pairs(groupby(vcf_files_df, :SAMPLES))
            put!(channel, (key.SAMPLES, group))
        end
    end
end

function split_bed_file_to_vcf(prefix; output_dir=".", samples_per_file=5_000)
    sample_batches = SequentialGWAS.make_sample_batches(prefix; output_dir=output_dir, samples_per_file=samples_per_file)
    chromosomes = SequentialGWAS.get_chromosomes(prefix)
    Threads.@threads for chr in chromosomes
        for samples_batch in sample_batches
            samples_ext = splitext(basename(samples_batch))[1]
            output_prefix = joinpath(output_dir, string(chr, "_", samples_ext))
            run(Cmd([
                "plink2",
                "--bfile", prefix,
                "--threads", "1",
                "--output-chr", "chr26",
                "--keep", samples_batch,
                "--chr", string(chr),
                "--export", "vcf-4.2", "id-delim=@",
                "--out", output_prefix
            ]))
            run(Cmd([
                "bgzip",
                "--force",
                string(output_prefix, ".vcf")
            ]))
        end
    end
end

function send_job_to_topmed(group, jobname, token, password; r2=0.8)
    cmd = Cmd([
        "curl", 
        "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver",
        "-X", "POST",
        "-H", "X-Auth-Token: $token",
        "-F", "job-name=$jobname",
        [string("-F files=@$(realpath(f))") for f in group.FILE]...,
        "-F", "refpanel=apps@topmed-r3",
        "-F", "build=hg38",
        "-F", "phasing=eagle",
        "-F", "password=$password",
        "-F", "population=all",
        "-F", "meta=yes",
        "-F", "r2Filter=$r2",
    ])
    job_details = JSON.parse(read(cmd, String))
    job_details["success"] == true || throw(error(job_details["message"]))
    return job_details
end

function get_job_status(job_id, token)
    cmd = Cmd([
        "curl", 
        "-H", "X-Auth-Token: $token", 
        string("https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/$job_id")
    ])
    return JSON.parse(read(cmd, String))
end

function wait_for_completion(token, job_id; rate=60)
    while true
        status = get_job_status(job_id, token)
        state = status["state"]
        if state == 5 || state == 6
            throw(error("Job $job_id failed."))
        elseif state == 4
            return status
        else
            sleep(rate)
        end
    end
end

function get_download_list(status, job_output_dir)
    download_list = []
    for output in status["outputParams"]
        if output["name"] !== "logfile"
            mkpath(joinpath(job_output_dir, output["name"]))
            for file_dict in output["files"]
                push!(
                    download_list, 
                    Dict(
                        "output_subdir" => output["name"],
                        "hash" => file_dict["hash"],
                        "filename" => file_dict["name"],
                    )
                )
            end
        end
    end
    return download_list
end

function download_files(download_list, job_output_dir, token, password)
    Threads.@threads for file_dict in download_list
        file_hash = file_dict["hash"]
        file_name = file_dict["filename"]
        download_subdir =  joinpath(job_output_dir, file_dict["output_subdir"])
        output_file = joinpath(download_subdir, file_name)
        # Download File
        run(Cmd([
                "curl", "-sL", 
                string("https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/$file_hash/$file_name"),
                "-H", "X-Auth-Token: $token",
                "-o", output_file
        ]))
        # Maybe Unzip
        if endswith(output_file, "zip")
            run(Cmd(["unzip", "-P", password, output_file, "-d", download_subdir]))
        end
    end
end

function download_results(status, token, password; output_dir=".")
    jobname = status["name"]
    job_output_dir = joinpath(output_dir, string("imputed_", jobname))
    download_list = SequentialGWAS.get_download_list(status, job_output_dir)
    download_files(download_list, job_output_dir, token, password)
end

function send_to_topmed_and_download(channel, token, password; refresh_rate=120, r2=0.8, output_dir=".")
    for (jobname, group) in channel
        job_details = send_job_to_topmed(group, jobname, token, password;r2=r2)
        job_id = job_details["id"]
        status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
        SequentialGWAS.download_results(status, token, password; output_dir=output_dir)
    end
end

function download_from_job_ids(jobs_file, token; password="abcde", refresh_rate=360, output_dir=".")
    job_ids = readlines(jobs_file)
    for job_id in job_ids
        status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
        SequentialGWAS.download_results(status, token, password; output_dir=output_dir)
    end
    return 0
end

function impute(genotypes_prefix, token_file; 
    password="abcde", 
    max_concurrent_submissions=3,
    refresh_rate=120,
    r2=0.8,
    samples_per_file=10_000,
    jobs_file=nothing,
    output_dir="."
    )
    token = read(token_file, String)
    if jobs_file !== nothing
        download_from_job_ids(jobs_file, token; password=password, refresh_rate=refresh_rate, output_dir=output_dir)
    else
        # Split the bed file into smaller VCF files for each chromosome
        vcfs_dir = joinpath(output_dir, "vcfs")
        mkpath(vcfs_dir)
        SequentialGWAS.split_bed_file_to_vcf(genotypes_prefix; output_dir=vcfs_dir, samples_per_file=samples_per_file)
        # Group files into a channel for submission
        vcf_files_channel = SequentialGWAS.get_vcf_files_channel(vcfs_dir)
        # Send for submission and download, maximum 3 concurrent tasks running at once on topmed
        tasks = [
            Threads.@spawn send_to_topmed_and_download(
                vcf_files_channel, 
                token, 
                password; 
                refresh_rate=refresh_rate,
                r2=r2,
                output_dir=output_dir
            ) for _ in 1:max_concurrent_submissions
        ]

        for task in tasks
            wait(task)
        end
        
    end
    return 0
end

function merge_imputed_genotypes_by_chr(chr, prefix)
    merge_list = "work/b8/d8c46db6b2c3c5cc6431722136609b/merge_list.txt"
    sample_list = "work/b8/d8c46db6b2c3c5cc6431722136609b/samples.txt"
    genotype_files = readlines(merge_list)
    samples = readlines(sample_list)
    samples_end = parse.(Int, [split(s, "-")[end] for s in samples])
    sorted_genotypes_and_samples = sort([(gf, se) for (gf, se) in zip(genotype_files, samples_end)], by= x->x[2])
    sorted_vcf_files = getindex.(sorted_genotypes_and_samples, 1)
    # Merge
    merge_output = read(Cmd([
            "bcftools", 
            "merge",
            "--threads", string(nthreads()),
            "-o", "$chr.merged.vcf.gz",
            "-O", "z",
            sorted_vcf_files...
            ]), 
        String
    )
    # Convert to BGEN
    read(Cmd([
        "plink2", 
        "--vcf", "$chr.merged.vcf.gz", 
        "--export", "bgen-1.2", 
        "--out", "$chr.merged.bgen"
        ]), 
        String
    )
end