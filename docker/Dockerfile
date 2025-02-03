FROM julia:1.11-bullseye

RUN apt-get update && apt-get install -y wget procps clang libblas-dev liblapack-dev libopenblas-dev tabix

# Install QCTOOL

RUN wget -O qctool.tgz https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.5-CentOS_Linux7.9-x86_64.tgz 
RUN tar -xvzf qctool.tgz
RUN mv qctool_v2.2.5-Ce-x86_64/qctool /usr/local/bin

# Install conda

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

RUN bash Miniforge3.sh -b -p /opt/miniforge3

RUN /opt/miniforge3/bin/mamba init

# Install regenie

RUN /opt/miniforge3/bin/mamba create -y -n regenie_env -c conda-forge -c bioconda regenie==4.0

# Install plink

RUN /opt/miniforge3/bin/mamba create -y -n plink_env -c conda-forge -c bioconda plink

# Install plink2

RUN /opt/miniforge3/bin/mamba create -y -n plink2_env -c conda-forge -c bioconda plink2

# Install UCSC liftover

RUN /opt/miniforge3/bin/mamba create -y -n liftover_env -c bioconda ucsc-liftover python=2.7

# Add Project code, install dependencies and build sysimage

COPY . /opt/sequential-gwas

RUN julia --project=/opt/sequential-gwas -e 'using Pkg; Pkg.instantiate()'

# RUN julia --project=/opt/sequential-gwas /opt/sequential-gwas/src/sysimage.jl
