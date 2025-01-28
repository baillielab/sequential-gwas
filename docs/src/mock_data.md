# Mock Data

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