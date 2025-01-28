# Development

This page contains some tips and tricks to help development. Because of the difficulty associated with development on ODAP it is assumed that most of the development will take place outside of it. For that purpose the original data has been mocked and is stored in `test/assets`. More information on the process can be found in [Mock Data](@ref)

## Build Docker Image

The image is built automatically during the continuous integration process and published on [Docker HUB](https://hub.docker.com/repository/docker/olivierlabayle/sequential-gwas/tags).

You can also build it locally if you have docker installed by running the following:

```bash
docker build -t sequential-gwas -f docker/Dockerfile .
```

(If running on MacOS with an ARM platform, add: `--platform linux/amd64`)

## Code Development on ODAP2

!!! warning 
    When you really need to.

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