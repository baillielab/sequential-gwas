version 1.0 

workflow merge_ukb_genomicc { 
    input { 
        String name 
    }

    call write_greeting { 
        input: greet_name = name 
    }
}

task write_greeting { 
    input {
        String greet_name 
    }

    command { 
        echo 'Hello, ${greet_name}!' > out.txt
    }

    output {
        File outfile = "out.txt"
    }
}