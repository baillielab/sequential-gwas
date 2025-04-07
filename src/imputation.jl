function write_sample_batches(prefix; output_prefix="genomicc", samples_per_file=5_000)
    fam = SequentialGWAS.read_fam(string(prefix, ".fam"))
    return map(Iterators.partition(1:nrow(fam), samples_per_file)) do indices
        filename = string(output_prefix, ".samples_", indices[1], "_", indices[end], ".keep")
        CSV.write(filename, fam[indices, [:FID, :IID]], delim=' ', header=false)
    end
end

function write_chromosome_list(genotypes_prefix; output_prefix="genomicc")
    bim = SequentialGWAS.read_bim(string(genotypes_prefix, ".bim"))
    open(string(output_prefix, ".chromosomes.txt"), "w") do io
        for chr in unique(bim.CHR_CODE)
            println(io, chr)
        end
    end
end

function get_vcf_files_channel(genotypes_prefix)
    dir = dirname(genotypes_prefix)
    dir = dir == "" ? "." : dir
    vcf_files = filter(x -> occursin(genotypes_prefix, x), readdir(dir, join=true))
    samples = [split(basename(f), ".")[end-2] for f in vcf_files]
    chromosomes = [split(basename(f), ".")[2] for f in vcf_files]
    vcf_files_df = DataFrame(FILE=vcf_files, SAMPLES=samples, CHR=chromosomes)
    return Channel() do channel
        for (key, group) in pairs(groupby(vcf_files_df, :SAMPLES))
            put!(channel, (key.SAMPLES, group))
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

function get_download_list(status)
    download_list = []
    for output in status["outputParams"]
        if output["name"] !== "logfile"
            for file_dict in output["files"]
                push!(
                    download_list, 
                    Dict(
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

function _download_topmed_file(file_dict, token, jobname; refresh_rate=360)
    file_hash = file_dict["hash"]
    file_name = file_dict["filename"]
    output_file = string(jobname, ".", file_name)
    while true
        read(Cmd([
            "curl", "-sL", 
            string("https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/$file_hash/$file_name"),
            "-H", "X-Auth-Token: $token",
            "-o", output_file
        ]), String)
        if has_download_suceeded(output_file)
            return output_file
        else
            sleep(refresh_rate)
        end
    end
end

function send_to_topmed_and_write_job_id(channel, token, password; refresh_rate=120, r2=0.8, output_prefix="genomicc")
    for (jobname, group) in channel
        job_details = send_job_to_topmed(group, jobname, token, password;r2=r2)
        job_id = job_details["id"]
        status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
        write(string(output_prefix, ".", job_id, ".txt"), job_id)
    end
end

function get_token(token_file)
    token = read(token_file, String)
    return endswith(token, "\n") ? token[1:end-1] : token
end

function get_download_list_and_checksum(job_id, token_file; refresh_rate=360)
    token = get_token(token_file)
    status = wait_for_completion(token, job_id; rate=refresh_rate)
    jobname = status["name"]
    download_list = get_download_list(status)
    # Download md5 checksums
    md5_index = only(findall(x -> x["filename"] == "results.md5", download_list))
    md5_file_dict = popat!(download_list, md5_index)
    _download_topmed_file(md5_file_dict, token, jobname; refresh_rate=refresh_rate)
    # Write files to be downloaded
    for file_dict in download_list
        open(string(jobname, ".todownload.", file_dict["filename"]), "w") do io
            JSON.print(io, file_dict)
        end
    end

    return 0
end

function get_md5_dict(md5_file)
    md5_df = CSV.read(md5_file, DataFrame; header=["MD5", "FILENAME"])
    return Dict(zip(md5_df.FILENAME, md5_df.MD5))
end

function download_topmed_file(job_id, token_file, file_info; md5_file=nothing, refresh_rate=360)
    token = SequentialGWAS.get_token(token_file)
    status = SequentialGWAS.wait_for_completion(token, job_id; rate=refresh_rate)
    jobname = status["name"]
    file_dict = open(JSON.parse, file_info)
    output_file = SequentialGWAS._download_topmed_file(file_dict, token, jobname; refresh_rate=refresh_rate)
    
    # Check hash for zip files
    if endswith(file_dict["filename"], "zip")
        expected_hash = SequentialGWAS.get_md5_dict(md5_file)[file_dict["filename"]]
        hash_to_file = readchomp(`md5sum $output_file`)
        file_hash = first(split(hash_to_file, " "))
        file_hash == expected_hash ||
            throw("Downloaded file's md5 checksum does not match expected md5 checksum. ($output_file)")
    end

    return 0
end

function write_imputation_split_lists(genotypes_prefix; output_prefix="genomicc", samples_per_file=20_000)
    write_sample_batches(genotypes_prefix; output_prefix=output_prefix, samples_per_file=samples_per_file)
    write_chromosome_list(genotypes_prefix; output_prefix=output_prefix)
end

function impute(genotypes_prefix, token_file; 
    password="abcde", 
    max_concurrent_submissions=3,
    refresh_rate=120,
    r2=0.8,
    output_prefix="genomicc"
    )
    token = get_token(token_file)
    # Split the bed file into smaller VCF files for each chromosome
    # Group files into a channel for submission
    vcf_files_channel = SequentialGWAS.get_vcf_files_channel(genotypes_prefix)
    # Send for submission and download, maximum 3 concurrent tasks running at once on topmed
    tasks = [
        Threads.@spawn send_to_topmed_and_write_job_id(
            vcf_files_channel, 
            token, 
            password; 
            refresh_rate=refresh_rate,
            r2=r2,
            output_prefix=output_prefix
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