def get_prefix(file){
    return file.toString().take(file.toString().lastIndexOf('.'))
}

def get_julia_cmd(cpus){
    def sysimageFile = new File("/opt/sequential-gwas/FlowOMMIC.so")
    if (workflow.profile == "dev") {
        return "julia --project=/opt/sequential-gwas --startup-file=no --threads=${cpus} /opt/sequential-gwas/bin/seq-gwas.jl"
    }
    else if (workflow.profile == "devsingularity"){
        return "JULIA_CPU_TARGET=generic JULIA_DEPOT_PATH=/tmp:\$JULIA_DEPOT_PATH julia --project=/opt/sequential-gwas --startup-file=no --threads=${cpus} /opt/sequential-gwas/bin/seq-gwas.jl"
    }
    else {
        return "JULIA_DEPOT_PATH=/tmp:\$JULIA_DEPOT_PATH julia --project=/opt/sequential-gwas --startup-file=no --threads=${cpus} --sysimage=${sysimageFile} /opt/sequential-gwas/bin/seq-gwas.jl"
    }        
}