FROM olivierlabayle/genomicc-base:main

RUN apt-get clean \
    && apt-get update 

# Create Julia App
COPY . /opt/sequential-gwas

# RUN cd /opt/sequential-gwas \
#     && julia --project --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

RUN cd /opt/sequential-gwas \
    && julia --startup-file=no bin/make_app.jl \
    && chmod +x build/bin/seq-gwas \
    && ./build/bin/seq-gwas --help

ENV PATH=/opt/sequential-gwas/build/bin/:${PATH}

# Update permission on /opt
RUN chmod -R 777 /opt


