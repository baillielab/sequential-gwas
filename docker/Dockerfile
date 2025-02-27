FROM olivierlabayle/genomicc-base:main

RUN apt-get clean \
    && apt-get update 

# Create Julia App
COPY . /opt/sequential-gwas

# RUN cd /opt/sequential-gwas \
#     && julia --project --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

ENV JULIA_DEPOT_PATH=/tmp:${PATH}

RUN cd /opt/sequential-gwas \
    && julia --startup-file=no bin/make_app.jl \
    && chmod -R 777 build \
    && ./build/bin/seq-gwas --help

ENV PATH=/opt/sequential-gwas/build/bin/:${PATH}




