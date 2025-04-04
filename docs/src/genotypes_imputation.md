# Genotypes Imputation

This workflow sends the genotypes for imputation to [TOPMed](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home), it follows the principles provided in their [documentation](https://topmedimpute.readthedocs.io/en/latest/).

!!! note "Platform"
    At the moment, it is impossible to run this workflow from ODAP because the TOPMed servers have not been white listed. Instead you can run it locally or from [eddie](https://digitalresearchservices.ed.ac.uk/resources/eddie). It is also highly beneficial to have 1 thread per chromosome file for downloads since these can be done concurrently.

## Inputs

For thiw workflow to work, you will need:

- Genotypes: likely the output of the [Combining Datasets](@ref) workflow.
- A TOPMed token: see [this page](http://topmedimpute.readthedocs.io/en/latest/api/#authentication).
- job ids (optional): Only to resume a crashed job.

Please have a look at the wokflow parameters below for how to setup the run.

## Running The Workflow

If the previous steps have been completed successfully you can run:

```bash
./run.sh ImputeAndDownload
```

If for some reason, all jobs get submitted but the workflow crashes before the results are downloaded, you can resume with:

```bash
./run.sh DownloadImputed
```

If not all jobs were submitted, it is likely better to cancel the running jobs manually on TOPMed and try the `ImputeAndDownload` again.

## Outputs

All outputs are produced in `PUBLISH_DIR` (defaults to `results`), the main outputs of the workflow are:

- `chr_P.{bgen,bgi,sample`: A set of imputed genotypes, one for each chromosome.

## Workflow Parameters

This is the list of all the pipeline's parameters, they can be set in the `run.config` file under the `params` section.

### Inputs Parameters

- `GENOTYPES_PREFIX`: Prefix to genotypes files.
- `TOPMED_TOKEN_FILE`: Path to the file containing your TOPMed API token.

### Important Options

- `TOPMED_ENCRYPTION_PASSWORD`: An encryption password
- `TOPMED_JOBS_LIST`: If the workflow crashes and you want to resume, list the job-ids in this file (one per line). Job ids can be obtained from the job url in TOPMed.
- `N_SAMPLES_PER_IMPUTATION_JOBS` (default: 10000): We can only send file of less than 200000 samples to TOPMed and the server only allows 3 jobs at a time. This number ideally splits your data in 3 roughly equal batches.
- `IMPUTATION_R2_FILTER` (default: 0.8): Only imputed variants passing the threshold are kept, set to 0 if you want to keep them all.

### Secondary Options

- `TOPMED_REFRESH_RATE` (default: 180): The frequency (in seconds) with which the workflow will monitor the imputation process to send further jobs.
- `TOPMED_MAX_PARALLEL_JOBS` (default: 3): The maximum number of concurrent imputation processes, this is limited to 3 at the moment by TOPMed
