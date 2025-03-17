# Prerequisites

## Platform: ODAP

Since the data is sensitive, it is meant to run on the [ODAP](https://odap.ac.uk/) platform (at least for now). To get access, you will need to contact the relevant person, at the moment [dominique.mccormick@ed.ac.uk](mailto:dominique.mccormick@ed.ac.uk). More information on accessing ODAP can be found [here](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide).

## Software

In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:

- Nextflow 24.10.3 ([see how to setup](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/nexflow))
- Singularity 3.9.4 (Should be ready to use)

## Importing Into ODAP

In order to upload the current state of the code, or a specific version of FlowOMICC into ODAP, we currently proceed via a shared folder with Dominique.

### Locally

The following assumes a specific git tag corresponding to a release for which a matching docker image exists, but the steps can be adapted to any need.

```bash
export FLOWOMICC_TAG="sha-9fd03f5"
```

In the shared folder.

1. Clone the repository or the relevant commit or tag:

```bash
git clone git@github.com:baillielab/sequential-gwas.git
git checkout $FLOWOMICC_TAG
```

2. Download and save the docker image
   
```bash
docker pull olivierlabayle/genomicc:$FLOWOMICC_TAG
docker save olivierlabayle/genomicc:$FLOWOMICC_TAG | gzip > genomicc.tar.gz
```

Then ask Dominique to uplaod the folder to ODAP.

### On ODAP

In the uploaded repository, build the singularity image:

```bash
singularity build genomicc.sif docker-archive:genomicc.tar.gz
```




