# ODAP

The original GenOMICC data is kept on the [ODAP](https://odap.ac.uk/) platform. To get access, you will need to contact the relevant person, at the moment [dominique.mccormick@ed.ac.uk](mailto:dominique.mccormick@ed.ac.uk). More information on accessing ODAP can be found [here](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide).

## Software

In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:

- Nextflow ([see how to setup](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/nexflow))
- Singularity (Should be ready to use)

## Importing Into ODAP

In order to upload the current state of the code, or a specific version of `genomicc-workflows` into ODAP, we currently proceed via a shared folder with Dominique.

### Locally

The following assumes a specific git tag corresponding to a release for which a matching docker image exists, but the steps can be adapted to any need.

```bash
export genomicc_workflows_tag="0.1.0"
```

In the shared folder.

1. Clone the repository or the relevant commit or tag:

```bash
git clone git@github.com:baillielab/genomicc-workflows.git
cd genomicc-workflows
git checkout v$genomicc_workflows_tag
```

2. Download and save the docker image
   
```bash
docker pull --platform linux/amd64 olivierlabayle/genomicc:$genomicc_workflows_tag
docker save olivierlabayle/genomicc:$genomicc_workflows_tag | gzip > genomicc.tar.gz
```

Then ask Dominique to uplaod the folder to ODAP.

### On ODAP

In the uploaded repository, build the singularity image:

```bash
gunzip genomicc.tar.gz && singularity build genomicc.sif docker-archive:genomicc.tar
```




