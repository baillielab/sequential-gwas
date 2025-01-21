FROM julia:1.11-bullseye

RUN apt-get update && apt-get install -y wget procps

# Install conda

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

RUN bash Miniforge3.sh -b -p /opt/miniforge3

RUN /opt/miniforge3/bin/mamba init

# Install regenie

RUN /opt/miniforge3/bin/mamba create -y -n regenie_env -c conda-forge -c bioconda regenie==4.0

# Install plink2

RUN /opt/miniforge3/bin/mamba create -y -n plink2_env -c conda-forge -c bioconda plink2

# Install UCSC liftover

RUN /opt/miniforge3/bin/mamba create -y -n liftover_env -c bioconda ucsc-liftover

# Add Project code and install dependencies

COPY . /opt/sequential-gwas

RUN julia --project=/opt/sequential-gwas -e 'using Pkg; Pkg.instantiate()'