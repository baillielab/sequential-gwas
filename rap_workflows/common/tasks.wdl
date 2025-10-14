version 1.0

import "structs.wdl"

task get_julia_cmd {
    input {
        String use_sysimage = "true"
        String threads = "auto"
    }
    command <<<
        julia_cmd_string="julia --project=/opt/genomicc-workflows --startup-file=no"
        if [[ "~{use_sysimage}" == "true" ]]; then
            julia_cmd_string+=" --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so"
        fi
        if [[ "~{threads}" == "auto" ]]; then
            julia_cmd_string+=" --threads=auto"
        fi
        julia_cmd_string+=" /opt/genomicc-workflows/bin/genomicc.jl"
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
    }
}

task ld_prune {
    input {
        String docker_image
        File high_ld_regions
        String chr
        File bed_file
        File bim_file
        File fam_file
        String output_prefix = "ld_pruned"
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
    }

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        plink2 \
            --bfile ${bed_prefix} \
            --indep-pairwise ~{ip_values}
        
        plink2 \
            --bfile ${bed_prefix} \
            --extract plink2.prune.in \
            --maf ~{maf} \
            --make-bed \
            --exclude range ~{high_ld_regions} \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "all",
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}