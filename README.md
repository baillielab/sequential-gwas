# sequential-gwas

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://baillielab.github.io/sequential-gwas/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://baillielab.github.io/sequential-gwas/dev/)
[![Build Status](https://github.com/baillielab/sequential-gwas/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/baillielab/sequential-gwas/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/baillielab/sequential-gwas/branch/main/graph/badge.svg)](https://codecov.io/gh/baillielab/sequential-gwas)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

The purpose of this project is to provide an end-to-end reproducible GWAS pipeline for the GenOMICC project. 

## ODAP

Since the data is sensitive, it is meant to run on the [ODAP](https://odap.ac.uk/) platform (at least for now). To get access, you will need to contact the relevant person, at the moment [dominique.mccormick@ed.ac.uk](dominique.mccormick@ed.ac.uk). More information on accessing ODAP can be found [here](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide).

## High Level Workflow Preview

TBD

## Input Data

Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This pipeline aims to make the most of all these diverse data sources.

### Genetic Data

#### Genotyping Arrays

| Brief Description | Period Begin | Period End | Genotyping Array | Genome Build | Directory | Genotypes File |
| :--------: | :------------: | :----------: | :----------------: | :------------: | :------------: | :------------: |
| Prehistoric r8 release | 04/05/2020   | 30/08/2021 | GSA-MD-24v3-0_A1               | GRC37        | wp5-gwas-r8-under90excl_2021Sep16 | PLINK_190921_0906/wp5-gwas-r8-under90excl_2021Sep16.ped |
| Before 2024 microarray | 20/09/2021   | 06/12/2023 | GSA-MD-24v3-0_A1 | GRC37 | 20210920_20231206_QC_VFinal | PLINK_040724_0954/20210920_20231206_QC_VFinal.ped |
| Since 2024 microarray  | 04/06/2024   | 10/06/2024 | GSA-MD-48v4-0_A1 | GRC38 | 20240604_20240610_QC_VFinal | PLINK_040724_0114/20240604_20240610_QC_VFinal.ped |

#### Whole Genome Sequencing

Currently the data available to me is a subset of ~500 000 SNPs split in files of groups of individuals in `/odp-beegfs/a015/linked_data/preqc/gws-genotyped`.

#### Erola's Data (Historical, TB Deleted)

As currently understood and relative to the data filesystem in `/odp-beegfs/a015/postqc/array/`:

- `imputation-output/chr*.info.gz` are outputs from TopMed imputation services from `erola_scripts/downloadTOPMED.sh`. These originate from merged/liftover genotyped files that have been lost.
- `imputation-output/Genomicc_chr*.txt` are probable outputs from `erola_scripts/filtervcf.sh` based on the above even though there seems to be an intermediate vcf file that was deleted.
- In `bgen/` are the likely output of `erola_scripts/vcftobgen.sh`.

### Phenotypes

## Build Docker Image Locally

```bash
docker build -t sequential-gwas -f docker/Dockerfile .
```

(If running on MacOS with an ARM platform, add: `--platform linux/amd64`)

## Dependencies

In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:

- Nextflow 24.10.3 ([see how to setup](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/nexflow))
- Singularity 3.9.4 (Should be ready to use)

## Usage

First import the project within ODAP, and, from the project root directory run:


```bash
nextflow run main.nf -profile odap
```

## For Developers: Code Development on ODAP2

First import the [docker image](https://hub.docker.com/repository/docker/olivierlabayle/sequential-gwas/general) as a singularity container within ODAP which we assume is called `sequential-gwas.sif`.

### Julia REPL

To get a shell while mounting the repo within the container:

```bash
singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif
```

```bash
singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif 
```

Then run the julia REPL

```bash
JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/mnt/sequential-gwas
```

### Extra Tools

Most tools are available within their conda environment, for instance regenie:

```bash
docker run -it --rm sequential-gwas /opt/miniforge3/bin/mamba run -n regenie_env regenie --help
```

(If running on MacOS with arm platform, add: `--platform linux/amd64`)
