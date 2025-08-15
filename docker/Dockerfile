FROM olivierlabayle/genomicc-base:scope

RUN apt-get clean \
    && apt-get update 

# Create Julia Sysimage

COPY . /opt/genomicc-workflows

RUN cd /opt/genomicc-workflows \
    && julia --project --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

RUN cd /opt/genomicc-workflows \
    && julia --project --startup-file=no --threads=auto bin/make_sysimage.jl



