def get_prefix(file){
    return file.toString().take(file.toString().lastIndexOf('.'))
}

def get_julia_cmd(cpus){
    def sysimageFile = new File("/opt/sequential-gwas/FlowOMMIC.so")
    def opts = workflow.profile in ["dev", "devsingularity"] ? "--threads=${cpus}" : "--sysimage=${sysimageFile} --threads=${cpus}"
    return "JULIA_DEPOT_PATH=/tmp:\$JULIA_DEPOT_PATH julia --project=/opt/sequential-gwas --startup-file=no ${opts} /opt/sequential-gwas/bin/seq-gwas.jl"
}