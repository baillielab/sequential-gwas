# Further Information

(To be moved/updated to an appropriate place)

The Illumina manifest files corresponding to each array can be downloaded from [Illumina's website](https://emea.support.illumina.com/array/array_kits/infinium-global-screening-array/downloads.html). The description of the manifest columns can be found [here](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565). A comparison of both the GSA-48v4 and GSA-24v3 arrays for the GRC38 genome build yields the following table (function `GenomiccWorkflows.OneTimeChecks.array_overlap`): 

| Description | Number of SNPs |
| :---: | :---: |
| In GSA-48v4 | 650 321 |
| In GSA-24v3 | 654 027 |
| Union | 702 515 |
| Intersection | 601 833 |
| Only in GSA-48v4 | 48 488 |
| Only in GSA-24v3 | 52 194 |

The arrays are thus quite similar with around 7% difference between the two and none being a subset of the other. This motivates the strategy to take the intersection of genotyping arrays before merging. Note that in reality the intersection might be smaller due to QC filtering of the input SNPs.

Furthermore, because the manifest files are too big to be version controlled, variants on the - RefStrand were extracted (function `GenomiccWorkflows.OneTimeChecks.make_snps_to_flip_list`) and stored as follows:
- GSA-MD-24v3-0_A1 (genome build GRC37) : `assets/GSA-24v3-0_A1-minus-strand.txt`
- GSA-MD-48v4-0_A1 (genome build GRC38) : `assets/GSA-48v4-0_20085471_D2-minus-strand.txt`

This has been done using the following command:

```bash
julia --project --startup-file=no bin/genomicc.jl snps-to-flip --help
```