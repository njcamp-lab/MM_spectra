# CoMMpass MM spectra derivation and modeling 

## Intro

Code used to generate spectra using the MMRF  CoMMpass multiple myeloma (MM) RNAseq data.  

more on the MMRF CoMMpass project : https://themmrf.org/finding-a-cure/our-work/the-mmrf-commpass-study/ 

***NOTE:*** files that end with the suffix `.nb.html` are R html notebook files that can be downloaded and opened in a web browser (they do not preview in GitHub) and contain output as well as code. The code is easily extracted from these files by clicking on the **Code** dropdown in the upper right hand corner of the notebook and selecting the *Download Rmd* option. 

## Input data 

The current CoMMpass data is in the NCI GDC data portal but is controlled access through dbGaP. 

https://portal.gdc.cancer.gov/projects/MMRF-COMMPASS 

https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000748.v7.p4

We used versions IA14 of the RNAseq transcript counts, patient data, and sample QC data that were previously available from the MMRF. This dataset is controlled access and we cannot include the dataset in this repo, see the links above for access.  

***NOTE:*** there are several file paths to input files that will have to be changed to reflect the file structure and location of the input files on your system should you choose to gain access and download them. The `vroom` and `read_csv` functions (from the `vroom` and `readr` R packages, respectively) are used for file input and require file paths. 

### Additional reference data 

We also used a version of `Homo_sapiens.GRCh37.74.gtf` that came with the CoMMpass data and was converted to a flat text csv file for gene annotations. 

We removed genes identified by Arora et al. (https://doi.org/10.1038/s41598-020-59516-z) as discordant using the gene list in the supplemental information of that paper. 

## Spectra derivation 

### Exploratory data analysis 

#### Batch variable and covariates for ComBat

It is important to correctly identify and specify the single batch variable for the ComBat algorithm that we use for batch correction (using the `ComBat` function from the `sva` packjage: https://www.bioconductor.org/packages/release/bioc/html/sva.html). There is a `Batch` variable in the sequencing QC information that comes with the CoMMpass data that we used for batch correction. 

We also include covariates of interest that are associated with batch (e.g by regression or Fisher's exact test). In the case of the CoMMpass project, we included: baseline vs. follow-up, age, gender, and all survival data (indicators and time). 

***Important note on `NA`s:*** the variables used in batch correcting with ComBat cannot contain `NA` values so all  people in the dataset that contain an `NA` for any variable used in the ComBat correction (either as the batch variable or covariates) must be removed from the analysis.  If there are `NA` values in the dataset ComBat will not run. 

#### Clinical variables for APV

To evaluate the added predictive value (APV, described later) of the spectra we needed to select variables for a clinical variable model to calculate the APV for spectra and compare spectra to other gene expression signatures. The variables used for this step should be those used to evaluate clinical risk. 

### Pre-processing

`CoMMpass_transcript_preprocessing_auto_followincl_final_20210601.R` contains the pre-processing steps to generate the input files required for the spectra derivation script, which are:

-  gene expression data
- gene metadata
- sample metadata 

These steps are outlined below.  

#### Expression data

The expression data file needs to have transcripts/genes as rows and samples as columns, with the first columns providing IDs and length information (names must match exactly)... 

*If using transcript count data* (as was done with CoMMpass): 

- `gene_id`
- `transcript_id`
- `length`

*If using gene count data* (without transcript information):

- `gene_id`
- `length`

Then all subsequent columns must be sample IDs and those sample IDs must match the rows of the sample metadata file (below). 

#### Gene metadata

This file provides information to the script about what genes to use in the spectra derivation.  

It has three columns (names must match exactly):

- `gene_id`
- `gene_name`
- `gene_status`

The gene IDs must match the gene IDs in the expression data.  

The gene names are used to match the discordant genes identified by Arora et al. 

Gene status tells the spectra derivation script how to treat each gene in the expression data and has the following possible values:  

- `allow`
- `force`
- `report`

`allow` is used for genes that should be evaluated to be included in spectra derivation - they will be evaluated by the gene QC steps and either kept in the derivation if they pass QC or excluded if they do not. In the CoMMpass data, protein coding genes on the autosomes were marked `allow` since those are the genes that we wanted to be included in the spectra derivation as long as they passed QC. 

`force` is used to include genes in the spectra derivation even if they fail QC. This is not often used in our experience, however, and example of why you might want to use this option is to include a well established cancer driver gene that is expressed in tumor and not in normal tissue in the case where you were including both tumor and normal tissue in the analysis. Depending on the balance of tumor and normal samples in your dataset, this gene may fail QC (for example if you have roughly equal number of samples of each type) but you might want to include it anyway due to the ability of the gene to distinguish tumor from normal tissue. No genes were marked as `force` for the CoMMpass project. 

`report` is used to include genes in the normalization process whether they make it through QC or not, but only include them in the spectra derivation if they make it through QC. For example, the genes in the UAMS and SBUK gene expression signatures were marked as `report` in the CoMMpass project since we needed the expression values to calculate the signatures but we didn't want them to be included in the spectra derivation if they did not pass QC. 

Genes that do not fall into any of the categories above should be removed from the gene metadata file

#### Sample metadata

This file provides information on which samples contribute to the spectra derivation and which are only normalized and batch corrected by ComBat. 

It has two required columns (names must match exactly):

- `sample_id`
- `status`

`status` has the following possible values: 

- `use`
- `correct_only`

`use` marks the sample for inclusion in the spectra derivation if the sample passes QC while `correct_only` samples will be included in ComBat correction and have their expression values reported if they pass QC but not included in the spectra derivation. In the CoMMpass data, baseline samples were marked `use` while `follow-up` samples were marked `correct_only` since we wanted to analyze the follow-up samples but not include them in the spectra derivation. 

If there is only the two required columns, no ComBat correction will be performed. 

If there is a third column, it will be treated as the batch variable for ComBat correction and any columns after the third will be used as covariates in the ComBat correction. Which variables should be included as batch and as covariates for ComBat is addressed above in the **exploratory data analysis** section. 

### Spectra derivation script 

The R script that processes the input data and derives the spectra is: `tissue_spectra_process.R`

To call the script from the command line the format is:  

`Rscript tissue_spectra_process.R <expression file name> <gene metadata file name> <sample metadata file name>`

where names inside `< >` are whatever you name the expression, gene metadata and sample metadate files. 

The output from this script includes (file name for CoMMpass data): 

- log info printed to stdout (`./SpectraData/logfile.txt`)
- a scree plot (`./SpectraData/pca_scree_plot-2021-06-01.pdf`)
- ComBat diagnostic plots (`./SpectraData/Diag/combat_diag_plot-2021-06-01.pdf`)
- normalized expression values (not included)
- spectra summary info (`./SpectraData/pca_details-2021-06-01.csv`)
- the rotation matrix and gene centers (`./SpectraData/spectra_rotation_matrix-gene_centers-2021-06-01.csv.gz`)
- spectra variables including the ComBat covariates (`./SpectraData/spectra_scores-batch_variables-2021-06-01.csv `)

### Post-processing 

The spectra variables are centered by default and we divide each by its standard deviation to put each spectrum on the same scale. All subsequent modeling steps were performed using these centered and scaled spectra variables. 

The standard deviation for the spectra can be derived from the variance (sd=sqrt(variance)) which is contained in `./SpectraData/pca_details-2021-06-01.csv` or derived directly from the scores (e.g. `CoMMpass_follow-up_samples_processing_20210602.nb.html`).  

The scaled spectra variables for the baseline samples are in `./SpectraData/Spectra_Variables_CoMMpass_Baseline_Samples_20210601.csv`

### Spectra description 

The bar chart of the variance for the spectra and the spectra barcode plot (Figure 3) were generated using: `plot_spectra_variance_barcodes.html`. 

## Predictive modeling

For each of the survival outcomes, we built initial predictive models listed here, then used bootstrap internal validation to estimate the optimism in our modeling process and adjust the estimates and confidence intervals (see next section). These files evaluate the GMM (using the `mclust` R package) for each outcome as well as make Kaplan-Meier survival curves (Figures 4 and 6). 

**OS** (overall survival): `CoMMpass_GMM_KM_curves_OS_20210602.nb.html`

**PFS** (progression free survival): `CoMMpass_GMM_KM_curves_PFS_20210604.nb.html`

**TTF** (time to treatment failure): `CoMMpass_GMM_KM_curves_TTF_20210603.nb.html`

The table of CoxPH beta values (Table 1) was generated using: `CoMMpass_survival_beta_table_20210604.nb.html` 

#### Bootstrap internal validation 

To mitigate overfitting, we used bootstrap internal validation (see *FJ Harrell, Regression Modeling Strategies, Springer 2015*). Each survival outcome has a separate bootstrapping script and analysis notebook. 

**OS** (overall survival, Table 2): `CoMMpass_OS_CoxPH_Spectra_UAMS_SBUK_bootstrap_production_version_20210601.R`, `CoMMpass_bootstrap_optimism_corrected_location-shifted_OS_HR_20210602.nb.html`

**PFS** (progression free survival): `CoMMpass_PFS_CoxPH_Spectra_UAMS_SBUK_bootstrap_production_version_20210604.R`, `CoMMpass_bootstrap_optimism_corrected_location-shifted_PFS_HR_20210604.nb.html`

**TTF** (time to treatment failure): `CoMMpass_TTF_CoxPH_Spectra_UAMS_SBUK_bootstrap_production_version_20210603.R`, `CoMMpass_bootstrap_optimism_corrected_location-shifted_TTF_HR_20210603.nb.html`

## Descriptive modeling 

The code for our description of the associations of the spectra with other variables (Figure 7) is contained in the R markdown file: `descriptive_models.html`. 

## Added predictive value (APV)

We used APV to describe how much information the spectra could contribute above and beyond a clinical model and to compare spectra with established gene expression signatures (see *FJ Harrell, Regression Modeling Strategies, Springer 2015* and https://www.fharrell.com/post/addvalue/)

APV was analyzed for each survival outcome separately:

**OS** (overall survival, Table 3): `CoMMpass_OS_RISS_age_compare_UAMS_SBUK_sigp_spectra_score_APV_20210602.nb.html`

**PFS** (progression free survival): `CoMMpass_PFS_RISS_age_compare_UAMS_SBUK_sigp_spectra_score_APV_20210604.nb.html`

**TTF** (time to treatment failure): `CoMMpass_TTF_RISS_age_compare_UAMS_SBUK_sigp_spectra_score_APV_20210603.nb.html`

## Follow-up samples 

Follow-up samples from CoMMpass that were corrected together with the baseline samples had their spectra variables calculated, analyzed, and visualized (Figure 5): `followup_spectra_OS.html`.

`CoMMpass_follow-up_samples_processing_20210602.nb.html` contains additional information about the follow-up samples and the standard deviations used to scale them. The standard deviation values themselves are in `./SpectraData/commpass_baseline_spectra_sds_20210602.csv`

The scaled spectra variables for the follow-up samples are in `./SpectraData/Spectra_Variables_CoMMpass_Followup_Samples.csv`

