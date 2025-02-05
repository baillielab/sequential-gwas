var documenterSearchIndex = {"docs":
[{"location":"cohorts/#Cohorts","page":"Cohorts Description","title":"Cohorts","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"The dataset is composed of 4 different cohorts: GenOMICC, MILD, REACT and ISARIC.","category":"page"},{"location":"cohorts/#GenOMICC-(Genetics-Of-Mortality-In-Critical-Care)","page":"Cohorts Description","title":"GenOMICC (Genetics Of Mortality In Critical Care)","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Link to the study webpage: https://genomicc.org/.","category":"page"},{"location":"cohorts/#GenOMICC-UK","page":"Cohorts Description","title":"GenOMICC UK","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Inclusion criteria:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Intensive Care Unit patients\nUK (mostly England and some specific recruitment centers if I understand correctly)","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Additional Information:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"The study also seems to be recruiting controls by matching in some way according to this page.\nNurses select who they think might be eligible which could produce cryptic bias as well.","category":"page"},{"location":"cohorts/#GenOMICC-International","page":"Cohorts Description","title":"GenOMICC International","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"TBD","category":"page"},{"location":"cohorts/#ISARIC","page":"Cohorts Description","title":"ISARIC","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Link to the study webpage: https://isaric.org/.","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Inclusion criteria:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Our subset of individuals all had Covid-19","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Additional Information:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"The severity varies and a severity score is provided\nMany biological and contextual covariates exist but are sparse","category":"page"},{"location":"cohorts/#MILD","page":"Cohorts Description","title":"MILD","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Link to the study webpage: ?","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Inclusion criteria:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Covid-19\nVolunteers (It is not clear in what way?)\nNot hospitalised ","category":"page"},{"location":"cohorts/#REACT-(Real-time-Assessment-of-Community-Transmission)","page":"Cohorts Description","title":"REACT (Real-time Assessment of Community Transmission)","text":"","category":"section"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Link to the study webpage https://www.imperial.ac.uk/medicine/research-and-impact/groups/react-study/.","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Inclusion criteria:","category":"page"},{"location":"cohorts/","page":"Cohorts Description","title":"Cohorts Description","text":"Covid-19\nIt is unclear if they were hospitalised\nIt is unclear what the selection process is","category":"page"},{"location":"julia_fns/#Index-of-Julia-Functions","page":"Index Of Julia Functions","title":"Index of Julia Functions","text":"","category":"section"},{"location":"julia_fns/","page":"Index Of Julia Functions","title":"Index Of Julia Functions","text":"Modules = [SequentialGWAS, SequentialGWAS.OneTimeChecks]","category":"page"},{"location":"julia_fns/#SequentialGWAS.read_bim-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.read_bim","text":"read_bim(prefix)\n\nColumns Description from: https://www.cog-genomics.org/plink/1.9/formats#bim\n\nChromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name\nVariant identifier\nPosition in morgans or centimorgans (safe to use dummy value of '0')\nBase-pair coordinate (1-based; limited to 231-2)\nAllele 1 (corresponding to clear bits in .bed; usually minor)\nAllele 2 (corresponding to set bits in .bed; usually major)\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.read_fam-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.read_fam","text":"read_fam(prefix)\n\nColumns Description from: https://www.cog-genomics.org/plink/1.9/formats#fam\n\nFamily ID ('FID')\nWithin-family ID ('IID'; cannot be '0')\nWithin-family ID of father ('0' if father isn't in dataset)\nWithin-family ID of mother ('0' if mother isn't in dataset)\nSex code ('1' = male, '2' = female, '0' = unknown)\nPhenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.report_qc_effect-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.report_qc_effect","text":"report_qc_effect(input_prefix, output_prefix)\n\nReport the effect of the QC filtering process on SNPs and samples.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.array_overlap-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.array_overlap","text":"array_overlap(manifest_file1, manifest_file2)\n\nTakes two manifest files from Illumina genotyping arrays and returns the intersection of the SNP names.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.identify_snps_to_flip-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.identify_snps_to_flip","text":"identify_snps_to_flip(manifest_file)\n\nAccording to this link,  the RefStrand column in the manifest file corresponds to the standard designation for all eukaryotic organisms used by HapMap and 1000 Genomes Project. Variants with RefStrand equal to - need to be flipped to the + strand.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.make_snps_to_flip_list-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.make_snps_to_flip_list","text":"make_snps_to_flip_list(output, manifest_file)\n\nTakes a manifest file from an Illumina genotyping array and writes a list of SNPs  that are on the - strand and need to be flipped to the + strand.\n\n\n\n\n\n","category":"method"},{"location":"outputs/#Output-Data","page":"Outputs","title":"Output Data","text":"","category":"section"},{"location":"outputs/","page":"Outputs","title":"Outputs","text":"TBD","category":"page"},{"location":"mock_data/#Mock-Data","page":"Mock Data","title":"Mock Data","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"In order to ease testing and development, we generate mock data that closely ressembles the original data while preserving the privacy of individuals in the cohorts. This page explains how this is done.","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"First but foremost, for all individuals, a mock sample ID is created. Since there is no sensitive information about individuals, these cannot be identified as long as the genetic data is perturbed. We explain below how this is done.","category":"page"},{"location":"mock_data/#Genotyping-Arrays","page":"Mock Data","title":"Genotyping Arrays","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Only a subset of genetic variations are kept (e.g. 100 out of 500 000)\nFor each individual, each variant is resampled independently from the cohort's empirical distribution. The probability of this operation to have no effect on an individual is difficult to estimate since it depends on each variant's alleles frequencies. If all variants were different it would be (1/n_samples)^nvariants which for the lower values of nsamples=1000 and nvariants=110, this is lower than 10^-300.","category":"page"},{"location":"mock_data/#Whole-Genome-Sequencing","page":"Mock Data","title":"Whole Genome Sequencing","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"The GVCF mock data arising from whole genome sequencing is built from a very small intersection of variants (e.g., 100) common to all genotyping arrays. Individuals are thus unidentifiable.","category":"page"},{"location":"mock_data/#Covariates","page":"Mock Data","title":"Covariates","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Since the covariates are not sensitive, the newly created odap identifier is simply forwarded to covariates.","category":"page"},{"location":"mock_data/#How-to-Mock","page":"Mock Data","title":"How to Mock","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"To run, on ODAP, assuming:","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"The data output by Dominique is in /odp-beefgs/a015/linked_data/preqc/array-pre-imputation/ and mounted in the singularity container in /mnt/data\nThe repo is mounted in /mnt/sequential-gwas (This is not necessary anymore once the code is in the container, just need to point to /opt/sequential-gwas)","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"singularity shell --bind /odp-beefgs/a015/linked_data/preqc/array-pre-imputation/:/mnt/data PATH_TO_SINGULARITY_IMAGE","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Then run ","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/opt/sequential-gwas /opt/sequential-gwas/bin/seq-gwas.jl","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Some patterns currently in the data but which should be removed by Dominique in the future are:","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"One individual has been genotyped twice\nSome sample ids have non-standard encoding","category":"page"},{"location":"development/#Development","page":"Development","title":"Development","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"This page contains some tips and tricks to help development. Because of the difficulty associated with development on ODAP it is assumed that most of the development will take place outside of it. For that purpose the original data has been mocked and is stored in test/assets. More information on the process can be found in Mock Data","category":"page"},{"location":"development/#Build-Docker-Image","page":"Development","title":"Build Docker Image","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"The image is built automatically during the continuous integration process and published on Docker HUB.","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"You can also build it locally if you have docker installed by running the following:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"docker build -t sequential-gwas -f docker/Dockerfile .","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"(If running on MacOS with an ARM platform, add: --platform linux/amd64)","category":"page"},{"location":"development/#Code-Development-on-ODAP2","page":"Development","title":"Code Development on ODAP2","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"When you really need to.","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"First import the docker image as a singularity container within ODAP which we assume is called sequential-gwas.sif.","category":"page"},{"location":"development/#Julia-REPL","page":"Development","title":"Julia REPL","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"To get a shell while mounting the repo within the container:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif ","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"Then run the julia REPL","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/mnt/sequential-gwas","category":"page"},{"location":"development/#Extra-Tools","page":"Development","title":"Extra Tools","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"Most tools are available within their conda environment, for instance regenie:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"docker run -it --rm sequential-gwas /opt/miniforge3/bin/mamba run -n regenie_env regenie --help","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"(If running on MacOS with arm platform, add: --platform linux/amd64)","category":"page"},{"location":"run/#Run-the-Workflow","page":"Run","title":"Run the Workflow","text":"","category":"section"},{"location":"run/#Platform:-ODAP","page":"Run","title":"Platform: ODAP","text":"","category":"section"},{"location":"run/","page":"Run","title":"Run","text":"Since the data is sensitive, it is meant to run on the ODAP platform (at least for now). To get access, you will need to contact the relevant person, at the moment dominique.mccormick@ed.ac.uk. More information on accessing ODAP can be found here.","category":"page"},{"location":"run/#Data-Access","page":"Run","title":"Data Access","text":"","category":"section"},{"location":"run/","page":"Run","title":"Run","text":"Once access to ODAP is given you will also need access to the Input Data. This is also handled by Dominique at the moment.","category":"page"},{"location":"run/#Software-Dependencies","page":"Run","title":"Software Dependencies","text":"","category":"section"},{"location":"run/","page":"Run","title":"Run","text":"In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:","category":"page"},{"location":"run/","page":"Run","title":"Run","text":"Nextflow 24.10.3 (see how to setup)\nSingularity 3.9.4 (Should be ready to use)","category":"page"},{"location":"run/#Run","page":"Run","title":"Run","text":"","category":"section"},{"location":"run/","page":"Run","title":"Run","text":"If the previous steps have been completed successfully you can simply run:","category":"page"},{"location":"run/","page":"Run","title":"Run","text":"nextflow run main.nf -profile odap -resume","category":"page"},{"location":"input_data/#Input-Data","page":"Input Data","title":"Input Data","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This page aims at making it clear what are the inputs and outputs to the pipeline.","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"All file paths are currently reported with respect to my a015 project but this will be standardized in the future.","category":"page"},{"location":"input_data/#Genetic-Data","page":"Input Data","title":"Genetic Data","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Depending on the time period, individuals may have been genotyped, whole genome sequenced or both.","category":"page"},{"location":"input_data/#Genotyping-Arrays","page":"Input Data","title":"Genotyping Arrays","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Brief Description Period Begin Period End Genotyping Array Genome Build Directory Genotypes File\nPrehistoric r8 release 04/05/2020 30/08/2021 GSA-MD-24v3-0_A1 GRC37 wp5-gwas-r8-under90excl_2021Sep16 PLINK1909210906/wp5-gwas-r8-under90excl_2021Sep16.ped\nBefore 2024 microarray 20/09/2021 06/12/2023 GSA-MD-24v3-0_A1 GRC37 2021092020231206QC_VFinal PLINK0407240954/2021092020231206QC_VFinal.ped\nSince 2024 microarray 04/06/2024 10/06/2024 GSA-MD-48v4-0_A1 GRC38 2024060420240610QC_VFinal PLINK0407240114/2024060420240610QC_VFinal.ped","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"The Illumina manifest files corresponding to each array can be downloaded from Illumina's website. The description of the manifest columns can be found here. A comparison of both the GSA-48v4 and GSA-24v3 arrays for the GRC38 genome build yields the following table (function SequentialGWAS.OneTimeChecks.array_overlap): ","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Description Number of SNPs\nIn GSA-48v4 650 321\nIn GSA-24v3 654 027\nUnion 702 515\nIntersection 601 833\nOnly in GSA-48v4 48 488\nOnly in GSA-24v3 52 194","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"The arrays are thus quite similar with around 7% difference between the two and none being a subset of the other. This motivates the strategy to take the intersection of genotyping arrays before merging. Note that in reality the intersection might be smaller due to QC filtering of the input SNPs.","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Furthermore, because the manifest files are too big to be version controlled, variants on the - RefStrand were extracted (function SequentialGWAS.OneTimeChecks.make_snps_to_flip_list) and stored as follows:","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"GSA-MD-24v3-0A1 (genome build GRC37) : `assets/GSA-24v3-0A1-minus-strand.txt`\nGSA-MD-48v4-0A1 (genome build GRC38) : `assets/GSA-48v4-020085471_D2-minus-strand.txt`","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"This has been done using the following command:","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"julia --project --startup-file=no bin/seq-gwas.jl snps-to-flip --help","category":"page"},{"location":"input_data/#Whole-Genome-Sequencing","page":"Input Data","title":"Whole Genome Sequencing","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"Currently the data available to me is a subset of ~500 000 SNPs split in files of groups of individuals in /odp-beegfs/a015/linked_data/preqc/gws-genotyped.","category":"page"},{"location":"input_data/#Erola's-Historical-Data-(This-paragraph-will-be-deleted)","page":"Input Data","title":"Erola's Historical Data (This paragraph will be deleted)","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"As currently understood and relative to the data filesystem in /odp-beegfs/a015/postqc/array/:","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"imputation-output/chr*.info.gz are outputs from TopMed imputation services from erola_scripts/downloadTOPMED.sh. These originate from merged/liftover genotyped files that have been lost.\nimputation-output/Genomicc_chr*.txt are probable outputs from erola_scripts/filtervcf.sh based on the above even though there seems to be an intermediate vcf file that was deleted.\nIn bgen/ are the likely output of erola_scripts/vcftobgen.sh.","category":"page"},{"location":"input_data/#Covariates-Data","page":"Input Data","title":"Covariates Data","text":"","category":"section"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"The covariate file is a CSV file in /odp-beegfs/a015/linked_data/postqc/array/a015_covariates.csv with the following columns:","category":"page"},{"location":"input_data/","page":"Input Data","title":"Input Data","text":"ODAP_ID: Unique invidual identifier (same as the IID/FID in the genetic data).\nCOHORTS: List of all cohorts the individual belongs to.\nAGEATRECRUITMENT: Age of the individual at recruitment. If an individual was recruited in two cohorts, typically (GenOMICC and ISARIC), GenOMICC takes precedence.\nREPORTED_SEX: Sex reported by an individual.\nPRIMARY_DIAGNOSIS: Only for GenOMICC patients and obtained via REDCap.\nDIAGNOSES: List of diagnoses inferred from REDCap\nISARICSEVERITYSCORE: Only for ISARIC individuals, severity of infection.\nGENETICSAMPLETYPE: Blood or Saliva.\nCALL_RATE: Genotyping call rate (saliva samples are believed to have lower quality).\nGENETICMEASUREMENTTECHNOLOGY: At this time either of (GSA-MD-24v3-0A1, GSA-MD-48v4-0A1, WGS).\nDATEATRECRUITMENT: Time marker of individual recruitement.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SequentialGWAS","category":"page"},{"location":"#SequentialGWAS","page":"Home","title":"SequentialGWAS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation for the SequentialGWAS pipeline within the Baillie lab, the repository is here.","category":"page"},{"location":"#Full-Workflow-DAG","page":"Home","title":"Full Workflow DAG","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: DAG)","category":"page"}]
}
