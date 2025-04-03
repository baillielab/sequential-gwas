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
            output_subdir = joinpath(job_output_dir, output["name"])
            mkpath(output_subdir)
            for file_dict in output["files"]
                push!(
                    download_list, 
                    Dict(
                        "output_subdir" => output_subdir,
                        "hash" => file_dict["hash"],
                        "filename" => file_dict["name"],
                    )
                )
            end
        end
    end
    return download_list
end

function has_download_suceeded(output_file)
    # If the number of max downloads are exceeded, then the content of the file will be a JSON with success=false
    first_line = readline(output_file)
    try
        content = JSON.parse(first_line)
        if content["success"] == false && content["message"] == "number of max downloads exceeded."
            @info "Could not download file $output_file due to max downloads exceeded. Waiting."
            return false
        end
    catch
        return true
    end
end

function download_file_from_topmed(file_dict, token, jobname; md5_dict = Dict(), refresh_rate=360)
    file_hash = file_dict["hash"]
    file_name = file_dict["filename"]
    download_subdir = file_dict["output_subdir"]
    output_prefix = string("samples-", jobname, ".")
    output_file = joinpath(download_subdir, string(output_prefix, file_name))
    while true
        read(Cmd([
            "curl", "-sL", 
            string("https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/$file_hash/$file_name"),
            "-H", "X-Auth-Token: $token",
            "-o", output_file
        ]), String)
        if has_download_suceeded(output_file)
            # Check hash is correct or redownload (This can happen because too many downloads happen at once and TOPMed does not like it)
            if file_name in keys(md5_dict)
                hash_to_file = readchomp(`md5sum $output_file`)
                file_md5 = first(split(hash_to_file, " "))
                if file_md5 !== md5_dict[file_name]
                    return output_file
                else
                    @info "Incorrect MD5 sum for file: $output_file, redownloading"
                    sleep(refresh_rate)
                end
            else
                return output_file
            end
        else
            sleep(refresh_rate)
        end
    end
end

function download_files(download_list, jobname, token, password; md5_dict=Dict(), refresh_rate=360)
    Threads.@threads for file_dict in download_list
        # Download File
        output_file = download_file_from_topmed(file_dict, token, jobname; md5_dict=md5_dict, refresh_rate=refresh_rate)
        # Maybe Unzip and rename with prefix
        if endswith(output_file, "zip")
            file_name = file_dict["filename"]
            download_subdir = file_dict["output_subdir"]
            output_prefix = string("samples-", jobname, ".")
            extract_dir = joinpath(download_subdir, replace(file_name, ".zip" => ""))
            run(Cmd(["unzip", "-P", password, output_file, "-d", extract_dir]))
            for extracted_file in readdir(extract_dir)
                mv(
                    joinpath(extract_dir, extracted_file), 
                    joinpath(download_subdir, string(output_prefix, extracted_file))
                )
            end
            rm(extract_dir)
        end
    end
end

function get_md5_map!(download_list, jobname, token, password; refresh_rate=refresh_rate)
    md5_index = only(findall(x -> x["filename"] == "results.md5", download_list))
    md5_file_dict = popat!(download_list, md5_index)
    output_file = download_file_from_topmed(md5_file_dict, token, jobname; refresh_rate=refresh_rate)
    md5_df = CSV.read(output_file, DataFrame; header=["MD5", "FILENAME"])
    return Dict(zip(md5_df.FILENAME, md5_df.MD5))
end

function download_results(status, token, password; output_dir=".", refresh_rate=360)
    jobname = status["name"]
    download_list = get_download_list(status, output_dir)
    md5_dict = get_md5_map!(download_list, jobname, token, password; refresh_rate=refresh_rate)
    download_files(download_list, jobname, token, password; md5_dict=md5_dict, refresh_rate=refresh_rate)
end

function send_to_topmed_and_write_job_id(channel, token, password; refresh_rate=120, r2=0.8, output_dir=".")
    for (jobname, group) in channel
        job_details = send_job_to_topmed(group, jobname, token, password;r2=r2)
        job_id = job_details["id"]
        status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
        write(joinpath(output_dir, string(job_id, ".txt")), job_id)
    end
end

function download_from_job_ids(jobs_file, token_file; password="abcde", refresh_rate=360, output_dir=".")
    token = read(token_file, String)
    job_ids = readlines(jobs_file)
    for job_id in job_ids
        status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
        SequentialGWAS.download_results(status, token, password; output_dir=output_dir, refresh_rate=refresh_rate)
    end
    return 0
end

function impute(genotypes_prefix, token_file; 
    password="abcde", 
    max_concurrent_submissions=3,
    refresh_rate=120,
    r2=0.8,
    samples_per_file=10_000,
    output_dir="."
    )
    token = read(token_file, String)
    # Split the bed file into smaller VCF files for each chromosome
    vcfs_dir = joinpath(output_dir, "vcfs")
    mkpath(vcfs_dir)
    SequentialGWAS.split_bed_file_to_vcf(genotypes_prefix; output_dir=vcfs_dir, samples_per_file=samples_per_file)
    # Group files into a channel for submission
    vcf_files_channel = SequentialGWAS.get_vcf_files_channel(vcfs_dir)
    # Send for submission and download, maximum 3 concurrent tasks running at once on topmed
    tasks = [
        Threads.@spawn send_to_topmed_and_write_job_id(
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

    return 0
end

function lower_sample_id(filename)
    sample_prefix = first(split(basename(filename), "."))
    return parse(Int, split(sample_prefix, "-")[end])
end

function merge_imputed_genotypes_by_chr(chr; output_prefix="$chr.merged")
    chr = "chr1"
    output_prefix="$chr.merged"
    merge_list = "merge_list.txt"
    vcf_files = readlines(merge_list)
    sort!(vcf_files, by=lower_sample_id)
    # Index
    @threads for vcf_file in vcf_files
        run(Cmd([
            "mamba", "run", "-n", "bcftools_env",
            "bcftools", "index", vcf_file
        ]))
    end
    # Merge
    merged_vcf = "$output_prefix.dose.vcf.gz"
    merge_output = read(Cmd([
            "mamba", "run", "-n", "bcftools_env",
            "bcftools", 
            "merge",
            "--threads", string(nthreads()),
            "-o", merged_vcf,
            "-O", "z",
            vcf_files...
            ]), 
        String
    )
    # Convert to BGEN
    read(Cmd([
        "plink2", 
        "--vcf", merged_vcf, 
        "--export", "bgen-1.2", 
        "--out", output_prefix
        ]), 
        String
    )
end