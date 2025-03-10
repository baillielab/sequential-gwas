def get_prefix(file){
    return file.toString().take(file.toString().lastIndexOf('.'))
}

def get_julia_cmd(cpus){
    return "julia --project=/opt/sequential-gwas --startup-file=no --threads ${cpus} /opt/sequential-gwas/bin/seq-gwas.jl"
}