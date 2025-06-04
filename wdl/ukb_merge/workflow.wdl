version 1.0 

import "liftover_genotypes.wdl" as liftover_genotypes

workflow merge_ukb_genomicc {
    input {
        Array[File]+ plink_beds
        Array[File]+ plink_bims
        Array[File]+ plink_fams
        File ucsc_chain
        File reference_fastagz
        String split_par_build_code = "hg38"
        Boolean output_autosomal = true
    }

    call liftover_genotypes.liftover_plink_beds { 
        input: hello_and_goodbye_input = "sub world" 
    }

    # call myTask { input: hello_and_goodbye.hello_output }

    output {
        String main_output = hello_and_goodbye.hello_output
    }
}

task filter_ukb_chr {
    input {
        File bgen_file
        File bgen_bgi_gile
        File bgen_sample_file
        File docker_image
        File samples_to_exclude
        File variants_to_include
    }

    command {
        qctool \
            -g ~{bgen_file} \
            -s ~{bgen_sample_file} \
            -excl-samples ~{samples_to_exclude} \
            -incl-variants ~{variants_to_include} \
            -og filtered.bgen
    }
}


task export_critical_table { 
    input {
        File dataset
        File fieldsfile
    }

    command {
        dx extract_dataset ~{dataset} \
        --fields-file ~{fieldsfile} \
        -o="critical_table.csv" \
        -icoding_option==RAW \
        -iheader_style=UKB-FORMAT \
        -ientity=hesin_critical \
        -ifield_names_file_txt=~{fieldnames}
    }

    runtime {
        dx_app: object {
            id: "applet-xxxx",
            type: "app" 
        }
        dx_timeout: "4H"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        File outfile = "critical.csv"
    }
}