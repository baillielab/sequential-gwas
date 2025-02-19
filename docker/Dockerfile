FROM julia:1.11-bullseye

ENV JULIA_DEPOT_PATH=/opt/julia-depot

ENV DEBIAN_FRONTEND=noninteractive

RUN rm /var/lib/dpkg/info/libc-bin.* \
    && apt-get clean \
    && apt-get update \
    && apt-get install -y \
        wget \
        procps \
        clang \
        libblas-dev \
        liblapack-dev \
        libopenblas-dev \
        tabix \
        unzip \
        openjdk-17-jdk \
        zlib1g-dev

# Install QCTOOL
RUN wget -O qctool.tgz https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.5-CentOS_Linux7.9-x86_64.tgz \
    && tar -xvzf qctool.tgz \
    && mv qctool_v2.2.5-Ce-x86_64/qctool /usr/local/bin

# KING 2.3.2
RUN wget https://www.kingrelatedness.com/KINGcode.tar.gz \
    && mkdir /opt/king \
    && tar -xzf KINGcode.tar.gz -C /opt/king \
    && cd /opt/king \
    && c++ -lm -lz -O2 -fopenmp -o king *.cpp \
    && cp king /usr/local/bin/

# Install GATK
RUN mkdir /opt/gatk \
    && wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip \
    && unzip gatk-4.6.1.0.zip -d /opt/

# Install conda

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3.sh -b -p /opt/miniforge3 \
    && /opt/miniforge3/bin/mamba init

# Install bcftools & samtools

RUN /opt/miniforge3/bin/mamba create -n bcftools_env -c conda-forge -c bioconda bcftools samtools

# Install regenie

RUN /opt/miniforge3/bin/mamba create -y -n regenie_env -c conda-forge -c bioconda regenie==4.0

# Install plink

RUN /opt/miniforge3/bin/mamba create -y -n plink_env -c conda-forge -c bioconda plink

# Install plink2

RUN /opt/miniforge3/bin/mamba create -y -n plink2_env -c conda-forge -c bioconda plink2

# Install UCSC liftover

RUN /opt/miniforge3/bin/mamba create -y -n liftover_env -c bioconda ucsc-liftover python=2.7

# Install Admixture

RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
    && tar -xvzf admixture_linux-1.3.0.tar.gz \
    && mv dist/admixture_linux-1.3.0/admixture /usr/local/bin

# Update Path

ENV PATH=/opt/miniforge3/bin/:/opt/gatk-4.6.1.0/:${PATH}

# Add Project code, install dependencies and build sysimage

COPY . /opt/sequential-gwas

RUN julia --project=/opt/sequential-gwas -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# RUN julia --project=/opt/sequential-gwas /opt/sequential-gwas/src/sysimage.jl

# RUN echo 'alias juliasys="julia --project=/opt/sequential-gwas -J/opt/sequential-gwas/SequentialGWAS.so"' >> ~/.bashrc
