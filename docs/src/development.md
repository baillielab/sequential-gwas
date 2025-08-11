# Development

This page contains some tips and tricks to help development. Because of the difficulty associated with development on ODAP it is assumed that most of the development will take place outside of it. For that purpose the original data has been mocked and is stored in `test/assets`. More information on the process can be found in [Mock Data](@ref)

## Build Docker Image

The image is built automatically during the continuous integration process and published on [Docker HUB](https://hub.docker.com/repository/docker/olivierlabayle/genomicc-workflows/tags).

You can also build it locally if you have docker installed by running the following:

```bash
docker build -t genomicc-workflows -f docker/Dockerfile .
```

(If running on MacOS with an ARM platform, add: `--platform linux/amd64`)

## Code Development on ODAP2

When you really need to.

First import the [docker image](https://hub.docker.com/repository/docker/olivierlabayle/genomicc-workflows/general) as a singularity container within ODAP which we assume is called `genomicc-workflows.sif`.

### Julia REPL

To get a shell while mounting the repo within the container:

```bash
singularity shell --bind $PWD:/mnt/genomicc-workflows genomicc-workflows.sif
```

```bash
singularity shell --bind $PWD:/mnt/genomicc-workflows genomicc-workflows.sif 
```

Then run the julia REPL

```bash
JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/mnt/genomicc-workflows
```

### Extra Tools

Most tools are available within their conda environment, for instance regenie:

```bash
docker run -it --rm genomicc-workflows /opt/miniforge3/bin/mamba run -n regenie_env regenie --help
```

(If running on MacOS with arm platform, add: `--platform linux/amd64`)


## Debugging on UKB RAP

### Cloud Workstation

To debug errors, it may be useful to run the code interactively, for this, you can use the [Cloud Workstation](https://documentation.dnanexus.com/developer/cloud-workstation). This [tutorial](https://academy.dnanexus.com/interactivecloudcomputing/cloudworkstation) may also be useful. To start an instance:

```bash
dx run \
  --instance-type mem2_ssd1_v2_x8 \
  -imax_session_length="10h" \
  -y \
  --ssh app-cloud_workstation
```

To import one or multiple files:

```bash
dx download file-J1P9y88JjZjXfq4Y5gYBxk86 file-J1P9y88JjZjbX18xFZ38qY2P
```

Then you can download the docker image and enter a container:

```bash
docker run -it --rm -v $PWD:/mnt/data olivierlabayle/genomicc:scope /bin/bash
```

The current directory is mounted to `/mnt/data`. From there, work as usual, for instance to start a Julia REPL:

```bash
julia --project=/opt/genomicc-workflows --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so --startup-file=no
```

Finally, when you are finished, terminate the job with the appropriate job-id:

```bash
dx terminate job-J1V4870JpYQP94jgb33y45qP
```