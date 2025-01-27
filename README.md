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

## Build Docker Image

The image is build automatically during the continuous integration process and published on [Docker HUB](https://hub.docker.com/repository/docker/olivierlabayle/sequential-gwas/tags).

You can build it locally if you have docker installed by running the following:

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

### Mock Data

In order to ease testing and development, we generate mock data that closely ressembles the original data while preserving the privacy of individuals in the cohoerts.

How this is done:

1. Sample IDs are all modified to a random odapxxxxx ID so that samples are unidentifiable.
2. Only a subset of genetic variations are kept (e.g. 100 out of 500 000)
3. For each individual, each variant is resampled independently from the cohort's empirical distribution. The probability of this operation to have no effect on an individual is difficult to estimate since it depends on each variant's alleles frequencies. If all variants were different it would be (1/n_samples)^nvariants which for the lower values of $nsamples=1000$ and $nvariants=110$, this is lower than $10^-300$.

To run, on ODAP, assuming:

- The data output by Dominique is in `/odp-beefgs/a015/linked_data/preqc/array-pre-imputation/` and mounted in the singularity container in `/mnt/data`
- The repo is mounted in `/mnt/sequential-gwas` (This is not necessary anymore once the code is in the container, just need to point to `/opt/sequential-gwas`)

```bash
singularity shell --bind /odp-beefgs/a015/linked_data/preqc/array-pre-imputation/:/mnt/data PATH_TO_SINGULARITY_IMAGE
```

Then run 
```bash
JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/opt/sequential-gwas /opt/sequential-gwas/bin/seq-gwas.jl
```

Some patterns currently in the data but which should be removed by Dominique in the future are:

- One individual has been genotyped twice
- Some sample ids have non-standard encoding