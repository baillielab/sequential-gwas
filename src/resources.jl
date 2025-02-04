function check_or_get_reference_genome(;
    path=joinpath("assets", "GRCh38.fa"), 
    url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
    verbosity=1)
    
    if !isfile(path)
        verbosity > 0 && @info(string("Downloading reference genome to: ", path))
        open(path, "w") do output_io
            HTTP.open("GET", url) do input_io
                for line in eachline(GzipDecompressorStream(BufferedInputStream(input_io)))
                    println(output_io, line)
                end
            end
        end
    else
        verbosity > 0 && @info(string("Reference genome already in : ", path))
    end
end

function check_or_get_illumina_GSA_24v3_GRC37_manifest(;
    path = joinpath("assets", "GSA_24v3_GRC37_manifest.csv"),
    url="https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A1-manifest-file-csv.zip",
    verbosity=1
    )
    if !isfile(path)
        verbosity > 0 && @info(string("Downloading GSA-24v3-0-A1-manifest-file-csv to: ", path))
        zip_file = Downloads.download(url)
        reader = ZipFile.Reader(zip_file)
        file = only(reader.files)
        write(path, read(file, String))
        close(reader)
    else
        verbosity > 0 && @info(string("GSA-24v3-0-A1-manifest-file-csv already in : ", path))
    end
end

function check_or_get_illumina_GSA_48v4_GRC38_manifest(;
    path = joinpath("assets", "GSA_48v4_GRC38_manifest.csv"),
    url="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/GSA-48v4-0_20085471_D2.csv",
    verbosity=1
    )
    if !isfile(path)
        verbosity > 0 && @info(string("Downloading GSA-48v4-0_20085471_D2.csv to: ", path))
        Downloads.download(url, path)
    else
        verbosity > 0 && @info(string("GSA-48v4-0_20085471_D2.csv already in : ", path))
    end
end