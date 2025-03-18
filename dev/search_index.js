var documenterSearchIndex = {"docs":
[{"location":"julia_fns/#Index-of-Julia-Functions","page":"Index Of Julia Functions","title":"Index of Julia Functions","text":"","category":"section"},{"location":"julia_fns/","page":"Index Of Julia Functions","title":"Index Of Julia Functions","text":"Modules = [SequentialGWAS, SequentialGWAS.OneTimeChecks]","category":"page"},{"location":"julia_fns/#SequentialGWAS.get_action-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.get_action","text":"Compares the variant to the 1000 GP information and returns an action to take together with a reason.\n\nPotential actions are\n\nDROP\nKEEP\nFLIP\n\nSome particularly unexpected ACTION (REASON) are:\n\n\"KEEP (REVERSE-REF-ALT)\"\n\"FLIP (COMPLEMENT-NOT-MATCHING-KGP)\"\n\nBecause they mean the minor/major alleles are reversed in our dataset as compared to the reference KGP.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.kgp_unrelated_individuals-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.kgp_unrelated_individuals","text":"Only keeps the first individual within each family and writes them to outfile.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.read_bim-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.read_bim","text":"read_bim(file)\n\nColumns Description from: https://www.cog-genomics.org/plink/1.9/formats#bim\n\nChromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name\nVariant identifier\nPosition in morgans or centimorgans (safe to use dummy value of '0')\nBase-pair coordinate (1-based; limited to 231-2)\nAllele 1 (corresponding to clear bits in .bed; usually minor)\nAllele 2 (corresponding to set bits in .bed; usually major)\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.read_fam-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.read_fam","text":"read_fam(file)\n\nColumns Description from: https://www.cog-genomics.org/plink/1.9/formats#fam\n\nFamily ID ('FID')\nWithin-family ID ('IID'; cannot be '0')\nWithin-family ID of father ('0' if father isn't in dataset)\nWithin-family ID of mother ('0' if mother isn't in dataset)\nSex code ('1' = male, '2' = female, '0' = unknown)\nPhenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.read_map-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.read_map","text":"read_map(file)\n\nColumns Description from: https://www.cog-genomics.org/plink/1.9/formats\n\nChromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.\nVariant identifier\nPosition in morgans or centimorgans (optional; also safe to use dummy value of '0')\nBase-pair coordinate\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.report_qc_effect-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.report_qc_effect","text":"report_qc_effect(input_prefix, output_prefix)\n\nReport the effect of the QC filtering process on SNPs and samples.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.write_map-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.write_map","text":"write_map(file_prefix, array)\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.write_release_samples_to_drop-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.write_release_samples_to_drop","text":"We drop duplicate individuals according to the following priority:\n\nWGS > More Recent Array > Older Array\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.array_overlap-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.array_overlap","text":"array_overlap(manifest_file1, manifest_file2)\n\nTakes two manifest files from Illumina genotyping arrays and returns the intersection of the SNP names.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.identify_snps_to_flip-Tuple{Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.identify_snps_to_flip","text":"identify_snps_to_flip(manifest_file)\n\nAccording to this link,  the RefStrand column in the manifest file corresponds to the standard designation for all eukaryotic organisms used by HapMap and 1000 Genomes Project. Variants with RefStrand equal to - need to be flipped to the + strand.\n\n\n\n\n\n","category":"method"},{"location":"julia_fns/#SequentialGWAS.OneTimeChecks.make_snps_to_flip_list-Tuple{Any, Any}","page":"Index Of Julia Functions","title":"SequentialGWAS.OneTimeChecks.make_snps_to_flip_list","text":"make_snps_to_flip_list(output, manifest_file)\n\nTakes a manifest file from an Illumina genotyping array and writes a list of SNPs  that are on the - strand and need to be flipped to the + strand.\n\n\n\n\n\n","category":"method"},{"location":"combining_with_ukb/#Combining-with-UK-Biobank","page":"Combining with UK Biobank","title":"Combining with UK Biobank","text":"","category":"section"},{"location":"misc/#Further-Information","page":"Further Information","title":"Further Information","text":"","category":"section"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"(To be moved/updated to an appropriate place)","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"The Illumina manifest files corresponding to each array can be downloaded from Illumina's website. The description of the manifest columns can be found here. A comparison of both the GSA-48v4 and GSA-24v3 arrays for the GRC38 genome build yields the following table (function SequentialGWAS.OneTimeChecks.array_overlap): ","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"Description Number of SNPs\nIn GSA-48v4 650 321\nIn GSA-24v3 654 027\nUnion 702 515\nIntersection 601 833\nOnly in GSA-48v4 48 488\nOnly in GSA-24v3 52 194","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"The arrays are thus quite similar with around 7% difference between the two and none being a subset of the other. This motivates the strategy to take the intersection of genotyping arrays before merging. Note that in reality the intersection might be smaller due to QC filtering of the input SNPs.","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"Furthermore, because the manifest files are too big to be version controlled, variants on the - RefStrand were extracted (function SequentialGWAS.OneTimeChecks.make_snps_to_flip_list) and stored as follows:","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"GSA-MD-24v3-0A1 (genome build GRC37) : `assets/GSA-24v3-0A1-minus-strand.txt`\nGSA-MD-48v4-0A1 (genome build GRC38) : `assets/GSA-48v4-020085471_D2-minus-strand.txt`","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"This has been done using the following command:","category":"page"},{"location":"misc/","page":"Further Information","title":"Further Information","text":"julia --project --startup-file=no bin/seq-gwas.jl snps-to-flip --help","category":"page"},{"location":"mock_data/#Mock-Data","page":"Mock Data","title":"Mock Data","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"In order to ease testing and development, we generate mock data that closely ressembles the original data while preserving the privacy of individuals in the cohorts. This page explains how this is done.","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"First but foremost, for all individuals, a mock sample ID is created. Since there is no sensitive information about individuals, these cannot be identified as long as the genetic data is perturbed. We explain below how this is done.","category":"page"},{"location":"mock_data/#Genotyping-Arrays","page":"Mock Data","title":"Genotyping Arrays","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Only a subset of genetic variations are kept (e.g. 100 out of 500 000)\nFor each individual, each variant is resampled independently from the cohort's empirical distribution. The probability of this operation to have no effect on an individual is difficult to estimate since it depends on each variant's alleles frequencies. If all variants were different it would be (1/n_samples)^nvariants which for the lower values of nsamples=1000 and nvariants=110, this is lower than 10^-300.","category":"page"},{"location":"mock_data/#Whole-Genome-Sequencing","page":"Mock Data","title":"Whole Genome Sequencing","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"The GVCF mock data arising from whole genome sequencing is built from a very small intersection of variants (e.g., 100) common to all genotyping arrays. Individuals are thus unidentifiable.","category":"page"},{"location":"mock_data/#Covariates","page":"Mock Data","title":"Covariates","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Since the covariates are not sensitive, the newly created odap identifier is simply forwarded to covariates.","category":"page"},{"location":"mock_data/#How-to-Mock","page":"Mock Data","title":"How to Mock","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"To run, on ODAP, assuming:","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"The data output by Dominique is in /odp-beefgs/a015/linked_data/preqc/array-pre-imputation/ and mounted in the singularity container in /mnt/data\nThe repo is mounted in /mnt/sequential-gwas (This is not necessary anymore once the code is in the container, just need to point to /opt/sequential-gwas)","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"singularity shell --bind /odp-beefgs/a015/linked_data/preqc/array-pre-imputation/:/mnt/data PATH_TO_SINGULARITY_IMAGE","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Then run ","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/opt/sequential-gwas /opt/sequential-gwas/bin/seq-gwas.jl","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"I also manually:","category":"page"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"Added duplicate sample IDs to reproduce what is in the data.\nChanged the position of GSA-rs114361133 in test/assets/genomicc/genotyping_arrays/mock.release_2021_2023.map to be unliftable.","category":"page"},{"location":"mock_data/#Thousands-Genomes","page":"Mock Data","title":"Thousands Genomes","text":"","category":"section"},{"location":"mock_data/","page":"Mock Data","title":"Mock Data","text":"The data was downloaded from the 1000GP FTP and pruned using the bin/make_thousand_genomes_filter_files.jl script.","category":"page"},{"location":"development/#Development","page":"Development","title":"Development","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"This page contains some tips and tricks to help development. Because of the difficulty associated with development on ODAP it is assumed that most of the development will take place outside of it. For that purpose the original data has been mocked and is stored in test/assets. More information on the process can be found in Mock Data","category":"page"},{"location":"development/#Build-Docker-Image","page":"Development","title":"Build Docker Image","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"The image is built automatically during the continuous integration process and published on Docker HUB.","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"You can also build it locally if you have docker installed by running the following:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"docker build -t sequential-gwas -f docker/Dockerfile .","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"(If running on MacOS with an ARM platform, add: --platform linux/amd64)","category":"page"},{"location":"development/#Code-Development-on-ODAP2","page":"Development","title":"Code Development on ODAP2","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"When you really need to.","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"First import the docker image as a singularity container within ODAP which we assume is called sequential-gwas.sif.","category":"page"},{"location":"development/#Julia-REPL","page":"Development","title":"Julia REPL","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"To get a shell while mounting the repo within the container:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"singularity shell --bind $PWD:/mnt/sequential-gwas sequential-gwas.sif ","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"Then run the julia REPL","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH:/root/.julia julia --project=/mnt/sequential-gwas","category":"page"},{"location":"development/#Extra-Tools","page":"Development","title":"Extra Tools","text":"","category":"section"},{"location":"development/","page":"Development","title":"Development","text":"Most tools are available within their conda environment, for instance regenie:","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"docker run -it --rm sequential-gwas /opt/miniforge3/bin/mamba run -n regenie_env regenie --help","category":"page"},{"location":"development/","page":"Development","title":"Development","text":"(If running on MacOS with arm platform, add: --platform linux/amd64)","category":"page"},{"location":"gwas/#GWAS","page":"GWAS","title":"GWAS","text":"","category":"section"},{"location":"gwas/#DAG","page":"GWAS","title":"DAG","text":"","category":"section"},{"location":"gwas/","page":"GWAS","title":"GWAS","text":"<iframe src=\"../assets/gwas_dag.html\" width=\"100%\" height=\"800px\"></iframe>","category":"page"},{"location":"genotypes_imputation/#Genotypes-Imputation","page":"Genotypes Imputation","title":"Genotypes Imputation","text":"","category":"section"},{"location":"combining_datasets/#Combining-Datasets","page":"Combining Datasets","title":"Combining Datasets","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"This workflow combines the various available data sources into a unified dataset.","category":"page"},{"location":"combining_datasets/#Inputs","page":"Combining Datasets","title":"Inputs","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"Since GenOMICC is an ongoing project where the data is continuously collected, it came and will continue to arrive in different formats. This page aims at making it clear what are the inputs and outputs to the pipeline.","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"All file paths are currently reported with respect to my a015 project but this will be standardized in the future.","category":"page"},{"location":"combining_datasets/#External-Resources","page":"Combining Datasets","title":"External Resources","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"As well as our in-house data, the pipeline depends on external reference data. In principle these files should already be present on ODAP and there is nothing you need to do.","category":"page"},{"location":"combining_datasets/#The-1000-GP","page":"Combining Datasets","title":"The 1000 GP","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"All VCF files and indexes present in this FTP folder\nThe associated 1000 GP pedigree file.","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"These files should be stored in the same folder which is defined by the KGP_DIR (default: /mnt/odap-beegfs/software/gwas-resources/1000-genomes-HC) Nextflow parameter.","category":"page"},{"location":"combining_datasets/#GATK","page":"Combining Datasets","title":"GATK","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The reference genome published by the Broad Institute.","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"This should be in a folder defined by the GATK_DIR (default: /mnt/odap-beegfs/software/gwas-resources/gatk) Nextflow parameter.","category":"page"},{"location":"combining_datasets/#Genetic-Data","page":"Combining Datasets","title":"Genetic Data","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The pipeline requires both genotyping arrays and whole-genome sequencing (wgs) data. These should be all be located within a pre-qc directory and organised as follows:","category":"page"},{"location":"combining_datasets/#Genotyping-Arrays","page":"Combining Datasets","title":"Genotyping Arrays","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"There are three filesets and three corresponding subfolders:","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The r8 release: Corresponds to genotyping data generated before 2021. The genotyping chip was the Illumina GSA-MD-24v3-0A1 and the associated genome build GRCh37. The corresponding subfolder is `wp5-gwas-r8-under90excl2021Sep16`.\nThe 2021-2023 release:  Corresponds to genotyping data generated between 2021 and 2023. It was also genotyped using the Illumina GSA-MD-24v3-0A1 chip and the genome build is also GRCh37. The corresponding subfolder is `2021092020231206QCVFinal`.\nThe 2024-now release: Corresponds to the latest fileset. The Illumina GSA-MD-48v4-0A1 chip was used and the genome build is GRCh38. The corresponding subfolder is `2024060420240610QCVFinal`.","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The following tables summarises the above","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"Brief Description Period Begin Period End Genotyping Array Genome Build Directory Genotypes File\nPrehistoric r8 release 04/05/2020 30/08/2021 GSA-MD-24v3-0_A1 GRC37 wp5-gwas-r8-under90excl_2021Sep16 PLINK1909210906/wp5-gwas-r8-under90excl_2021Sep16.ped\nBefore 2024 microarray 20/09/2021 06/12/2023 GSA-MD-24v3-0_A1 GRC37 2021092020231206QC_VFinal PLINK0407240954/2021092020231206QC_VFinal.ped\nSince 2024 microarray 04/06/2024 10/06/2024 GSA-MD-48v4-0_A1 GRC38 2024060420240610QC_VFinal PLINK0407240114/2024060420240610QC_VFinal.ped","category":"page"},{"location":"combining_datasets/#Whole-Genome-Sequencing","page":"Combining Datasets","title":"Whole Genome Sequencing","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The wgs GVCF files are all located in a wgs-reheadered folder.","category":"page"},{"location":"combining_datasets/#Covariates-Data","page":"Combining Datasets","title":"Covariates Data","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"The covariate file is a CSV file in the pre-qc directory and named covariates.csv with the following columns:","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"ODAP_ID: Unique invidual identifier (same as the IID/FID in the genetic data).\nCOHORTS: List of all cohorts the individual belongs to.\nAGEATRECRUITMENT: Age of the individual at recruitment. If an individual was recruited in two cohorts, typically (GenOMICC and ISARIC), GenOMICC takes precedence.\nREPORTED_SEX: Sex reported by an individual.\nPRIMARY_DIAGNOSIS: Only for GenOMICC patients and obtained via REDCap.\nDIAGNOSES: List of diagnoses inferred from REDCap\nISARICSEVERITYSCORE: Only for ISARIC individuals, severity of infection.\nGENETICSAMPLETYPE: Blood or Saliva.\nCALL_RATE: Genotyping call rate (saliva samples are believed to have lower quality).\nGENETICMEASUREMENTTECHNOLOGY: At this time either of (GSA-MD-24v3-0A1, GSA-MD-48v4-0A1, WGS).\nDATEATRECRUITMENT: Time marker of individual recruitement.","category":"page"},{"location":"combining_datasets/#Running-The-Workflow","page":"Combining Datasets","title":"Running The Workflow","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"If the previous steps have been completed successfully you can run:","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"bash run.sh","category":"page"},{"location":"combining_datasets/#Pipeline-parameters","page":"Combining Datasets","title":"Pipeline parameters","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"This is the list of all the pipeline's parameters. In principle they don't need to be changed if the conventions in this documentation have been respected and are up to date. Otherwise, please feel free to open an issue.","category":"page"},{"location":"combining_datasets/#Output-Directories-Parameters","page":"Combining Datasets","title":"Output Directories Parameters","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"PUBLISH_DIR (default: \"results\"): Top level directory where outputs will be output.\nKGP_PUBLISH_DIR (default: \"results/kgp\"): Where 1000 Genome Project data will be output.\nARRAY_GENOTYPES_PUBLISH_DIR (default: \"results/array_genotypes\"): Where data associated with the genotyping arrays will be output.\nWGS_PUBLISH_DIR (default: \"results/wgs\"): Where data associated with the whole-genome sequencing data will be output.\nGATK_PUBLISH_DIR (default: \"results/gatk\"): Where data associated with GATK requirements will be output.\nMERGED_PUBLISH_DIR (default: \"results/merged\"): Where the merged genetic data will be output.","category":"page"},{"location":"combining_datasets/#External-Inputs-Parameters","page":"Combining Datasets","title":"External Inputs Parameters","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"RESOURCES_DIR (default: ./assets/resources\"): Path to all external resources.\nKGP_DIR (default: \"./assets/kgp\"): Path to the 1000 Genome Project specific resources.\nGATK_DIR (default: \"./assets/gatk\"): Path to GATK specific resources.\nGRC37_TO_GRC38_CHAIN_FILE (default: \"./assets/hg19ToHg38.over.chain.gz\"): Path to chain file used to liftover the GRCh37 genotypes to GRCh38.","category":"page"},{"location":"combining_datasets/#Basic-QC-Parameters","page":"Combining Datasets","title":"Basic QC Parameters","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"QC_GENOTYPE_MISSING_RATE (default: 0.1): Maximum missing rate per variant across all individuals. Variants above the threshold are dropped.\nQC_INDIVIDUAL_MISSING_RATE (default: 0.1): Maximum missing rate per individual across genotypes. Individuals above the threshold are dropped.\nQC_HWE_P (default: 1e-5): Used to identify potential technical artifacts and drop variants.\nQC_HWE_K (default: 0.001): Used together with QC_HWE_P","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"## Current Limitations","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"These are current limitations of the aggregation workflow:","category":"page"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"Only chromosomes 1 to 22 are processed.\nOnly bi-allelic SNPs are used.","category":"page"},{"location":"combining_datasets/#DAG","page":"Combining Datasets","title":"DAG","text":"","category":"section"},{"location":"combining_datasets/","page":"Combining Datasets","title":"Combining Datasets","text":"<iframe src=\"../assets/combining_datasets_dag.html\" width=\"100%\" height=\"2000px\"></iframe>","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SequentialGWAS","category":"page"},{"location":"#FlowOMICC","page":"Home","title":"FlowOMICC","text":"","category":"section"},{"location":"#Purpose","page":"Home","title":"Purpose","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This repository provides a set of Nextflow pipelines that can be used to analyse data generated by the GenOMICC project.","category":"page"},{"location":"#Cohorts","page":"Home","title":"Cohorts","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Even though our main dataset consists of GenoMICC data, for historical reasons (COVID-19 pandemic), we have a total of 4 different cohorts: GenOMICC, MILD, REACT and ISARIC. MILD, REACT and ISARIC were integrated in order to provide non severe \"controls\" for various COVID-19 susceptibility studies.","category":"page"},{"location":"#GenOMICC-(Genetics-Of-Mortality-In-Critical-Care)","page":"Home","title":"GenOMICC (Genetics Of Mortality In Critical Care)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Link to the study webpage: https://genomicc.org/.","category":"page"},{"location":"#GenOMICC-UK","page":"Home","title":"GenOMICC UK","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Inclusion criteria:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Intensive Care Unit patients\nUK (mostly England and some specific recruitment centers if I understand correctly)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Additional Information:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The study also seems to be recruiting controls by matching in some way according to this page.\nNurses select who they think might be eligible which could produce cryptic bias as well.","category":"page"},{"location":"#GenOMICC-International","page":"Home","title":"GenOMICC International","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is still being deployed and is not part of our current data.","category":"page"},{"location":"#ISARIC","page":"Home","title":"ISARIC","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Link to the study webpage: https://isaric.org/.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Inclusion criteria:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Our subset of individuals all had Covid-19","category":"page"},{"location":"","page":"Home","title":"Home","text":"Additional Information:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The severity varies and a severity score is provided\nMany biological and contextual covariates exist but are sparse","category":"page"},{"location":"#MILD","page":"Home","title":"MILD","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Link to the study webpage: ?","category":"page"},{"location":"","page":"Home","title":"Home","text":"Inclusion criteria:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Covid-19\nVolunteers (It is not clear in what way?)\nNot hospitalised ","category":"page"},{"location":"#REACT-(Real-time-Assessment-of-Community-Transmission)","page":"Home","title":"REACT (Real-time Assessment of Community Transmission)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Link to the study webpage https://www.imperial.ac.uk/medicine/research-and-impact/groups/react-study/.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Inclusion criteria:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Covid-19\nIt is unclear if they were hospitalised\nIt is unclear what the selection process is","category":"page"},{"location":"#Available-Workflows","page":"Home","title":"Available Workflows","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"These workflows will typically be run sequentially.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Combining Datasets\nGenotypes Imputation\nCombining with UK Biobank\nGWAS","category":"page"},{"location":"prerequisites/#Prerequisites","page":"Prerequisites","title":"Prerequisites","text":"","category":"section"},{"location":"prerequisites/#Platform:-ODAP","page":"Prerequisites","title":"Platform: ODAP","text":"","category":"section"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"Since the data is sensitive, it is meant to run on the ODAP platform (at least for now). To get access, you will need to contact the relevant person, at the moment dominique.mccormick@ed.ac.uk. More information on accessing ODAP can be found here.","category":"page"},{"location":"prerequisites/#Software","page":"Prerequisites","title":"Software","text":"","category":"section"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"In order to run the workflows in this repository only 2 software need to be installed, and they should already be present on ODAP:","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"Nextflow 24.10.3 (see how to setup)\nSingularity 3.9.4 (Should be ready to use)","category":"page"},{"location":"prerequisites/#Importing-Into-ODAP","page":"Prerequisites","title":"Importing Into ODAP","text":"","category":"section"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"In order to upload the current state of the code, or a specific version of FlowOMICC into ODAP, we currently proceed via a shared folder with Dominique.","category":"page"},{"location":"prerequisites/#Locally","page":"Prerequisites","title":"Locally","text":"","category":"section"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"The following assumes a specific git tag corresponding to a release for which a matching docker image exists, but the steps can be adapted to any need.","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"export FLOWOMICC_TAG=\"sha-9fd03f5\"","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"In the shared folder.","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"Clone the repository or the relevant commit or tag:","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"git clone git@github.com:baillielab/sequential-gwas.git\ngit checkout $FLOWOMICC_TAG","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"Download and save the docker image","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"docker pull olivierlabayle/genomicc:$FLOWOMICC_TAG\ndocker save olivierlabayle/genomicc:$FLOWOMICC_TAG | gzip > genomicc.tar.gz","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"Then ask Dominique to uplaod the folder to ODAP.","category":"page"},{"location":"prerequisites/#On-ODAP","page":"Prerequisites","title":"On ODAP","text":"","category":"section"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"In the uploaded repository, build the singularity image:","category":"page"},{"location":"prerequisites/","page":"Prerequisites","title":"Prerequisites","text":"singularity build genomicc.sif docker-archive:genomicc.tar.gz","category":"page"}]
}
