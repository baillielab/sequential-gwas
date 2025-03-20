FROM olivierlabayle/genomicc-base:main

RUN apt-get clean \
    && apt-get update 

# Create Julia Sysimage

COPY . /opt/sequential-gwas

RUN cd /opt/sequential-gwas \
    && julia --project --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

RUN cd /opt/sequential-gwas \
    && julia --project --startup-file=no --threads=auto bin/make_sysimage.jl



