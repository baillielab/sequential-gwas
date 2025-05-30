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
