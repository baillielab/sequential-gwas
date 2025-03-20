def get_prefix(file){
    return file.toString().take(file.toString().lastIndexOf('.'))
}

def get_julia_cmd(cpus){
    def sysimageFile = new File("/opt/sequential-gwas/FlowOMMIC.so")
    def opts = sysimageFile.exists() ? "--sysimage=${sysimageFile} --threads=${cpus}" : "--threads=${cpus}"
    return "JULIA_DEPOT_PATH=/tmp:\$JULIA_DEPOT_PATH julia --project=/opt/sequential-gwas --startup-file=no ${opts} /opt/sequential-gwas/bin/seq-gwas.jl"
}