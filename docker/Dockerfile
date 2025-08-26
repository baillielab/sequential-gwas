FROM olivierlabayle/genomicc-base:scope

RUN apt-get clean \
    && apt-get update 

# Create Julia Sysimage

COPY . /opt/genomicc-workflows

RUN cd /opt/genomicc-workflows \
    && julia --project --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

RUN cd /opt/genomicc-workflows \
    && /opt/miniforge3/bin/conda run -n R_env R_HOME=/opt/miniforge3/envs/R_env/lib/R julia --project --startup-file=no --threads=auto bin/make_sysimage.jl



