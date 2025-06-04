# Working with DNANexus

Perhaps surprisingly, most of the interactions with the UK Biobank RAP is made  directly from your local machine using the dx-toolkit.

## Installing Dependencies

### DNA Nexus Toolkit

First, install the dx-toolkit on your local computer by following the instructions on [this page](https://documentation.dnanexus.com/downloads). You will also need to download the compiler.

For reference, a quickstart guide to the toolkit can also be found [here](https://documentation.dnanexus.com/getting-started/cli-quickstart).

### Installing miniwdl

If you want to check the workflow syntax and run the `make check` operations, you can also download [miniwdl](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl).

### Cromwell

To create apps you will need [cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) as well.


## Logging in

Then you need to login to your UKB RAP account using:

```bash
dx login
```

You will be prompted for your username and password.

## Preparing the assets

```bash
docker save olivierlabayle/genomicc:main | gzip > assets/rap/genomicc.tar.gz
```

## Uploading Assets

1. Create an `assets` folder on the UKB RAP

2. Bring the GenOMICC data in the assets folder

   - 

3. Run:

```bash
export PROJECT_ID=project-J0pkqyQJpYQ133JG1p2J1qzv
dx upload -r assets/rap/ --destination $PROJECT_ID:/assets/
```

3. Uploading the docker image

The Docker iamge seems to be too big for `dx upload` and we need to use the [upload agent](https://documentation.dnanexus.com/downloads#installing-the-upload-agent) (`ua`):

```bash
export AUTH_TOKEN=XXX
ua --project $PROJECT_ID$ --auth-token $AUTH_TOKEN$ --folder /assets assets/rap/genomicc.tar.gz
```

## Extract Datasets

