FROM julia:1.11-bullseye

RUN apt-get update && apt-get install -y wget

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

RUN bash Miniforge3.sh -b

RUN /root/miniforge3/bin/mamba init

RUN /root/miniforge3/bin/mamba create -y -n regenie_env -c conda-forge -c bioconda regenie==4.0
