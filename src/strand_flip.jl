function identify_assay_and_control_blocks(manifest_file)
    assay_line = 0
    for (line_id, line) in enumerate(eachline(manifest_file))
        if line == "[Assay]"
            assay_line = line_id
        end
        if line == "[Controls]"
            assay_line !== 0 || throw(ArgumentError("No [Assay] section found in manifest file: $manifest_file"))
            return assay_line, line_id
        end
    end
    throw(ArgumentError("No [Controls] section found in manifest file: $manifest_file"))
end

function load_illumina_manifest_file(manifest_file)
    assay_line, controls_line = identify_assay_and_control_blocks(manifest_file)
    return CSV.read(manifest_file, DataFrame; 
        header=assay_line+1, 
        delim=",",
        limit=controls_line-assay_line-2,
        ntasks=1)
end

"""

According to [this link](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565), 
the RefStrand column in the manifest file corresponds to the standard designation for all eukaryotic organisms used by HapMap and 1000 Genomes Project.
Variants with RefStrand equal to `-` need to be flipped to the + strand.
"""
function identify_snps_to_flip(manifest_file)
    manifest = load_illumina_manifest_file(manifest_file)
    return filter(:RefStrand => ==("-"), manifest).Name
end

function make_snps_to_flip_list(output, manifest_file)
    snps_to_flip = identify_snps_to_flip(manifest_file)
    open(output, "w") do io
        for snp in snps_to_flip
            println(io, snp)
        end
    end
end
