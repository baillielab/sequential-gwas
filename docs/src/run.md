# Run the Workflow

## Platform: ODAP

Since the data is sensitive, it is meant to run on the [ODAP](https://odap.ac.uk/) platform (at least for now). To get access, you will need to contact the relevant person, at the moment [dominique.mccormick@ed.ac.uk](mailto:dominique.mccormick@ed.ac.uk). More information on accessing ODAP can be found [here](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide).

## Data Access

Once access to ODAP is given you will also need access to the [Input Data](@ref). This is also handled by Dominique at the moment.

## Software Dependencies

In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:

- Nextflow 24.10.3 ([see how to setup](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/nexflow))
- Singularity 3.9.4 (Should be ready to use)

## Run

If the previous steps have been completed successfully you can simply run:

```bash
nextflow run main.nf -profile odap -resume
```

