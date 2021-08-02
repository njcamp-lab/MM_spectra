# CoMMpass MM spectra derivation and modeling 

## intro

Code used to generate spectra using the MMRF  CoMMpass multiple myeloma (MM) RNAseq data.  

more on the MMRF CoMMpass project : https://themmrf.org/finding-a-cure/our-work/the-mmrf-commpass-study/ 



## input data 

The current CoMMpass data is in the NCI GDC data portal but is controlled access through dbGaP. 

https://portal.gdc.cancer.gov/projects/MMRF-COMMPASS 

https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000748.v7.p4

We used versions IA14 of the RNAseq transcript counts, patient data, and sample QC data that were previously available from the MMRF. 

***NOTE:*** there are several file paths to input files that will have to be changed to reflect the file structure and location of the input files on your system. `vroom` and `read_csv` (from the `vroom` and `readr` packages respectively) are used for file input and require file paths. 

### additional reference data 

We also used a version of `Homo_sapiens.GRCh37.74.gtf` that came with the CoMMpass data and was converted to a flat text csv file (deposited here as `Homo_sapiens.GRCh37.74.gtf.flat.csv`) 

We removed genes identified by Arora et al. (https://doi.org/10.1038/s41598-020-59516-z) as discordant using the gene list in the supplemental information of that paper. 

## spectra derivation 

### exploratory data analysis 

#### batch variable and covariates for ComBat



#### clinical variables for APV



### pre-processing

`CoMMpass_transcript_preprocessing_auto_followincl_final_20210601.R` contains the pre-processing steps to generate the input files required for the spectra derivation script, which are:

-  gene expression data
- gene metadata
- sample metadata 

These steps are outlined below.  

#### expression data

The expression data file needs to have transcripts/genes as rows and samples as columns, with the first columns providing IDs and length information (names must match exactly)... 

*If using transcript count data* (as was done with CoMMpass): 

- `gene_id`
- `transcript_id`
- `length`

*If using gene count data* (without transcript information):

- `gene_id`
- `length`

Then all subsequent columns must be sample IDs and those sample IDs must match the rows of the sample metadata file (below). 

#### gene metadata

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

#### sample metadata

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

### spectra derivation script 

The R script that processes the input data and derives the spectra is: `tissue_spectra_process.R`

The output from this script includes (file name for CoMMpass data: 

- log info printed to stdout (`./SpectraData/logfile.txt`)
- a scree plot (`./SpectraData/pca_scree_plot-2021-06-01.pdf`)
- ComBat diagnostic plots (`./SpectraData/Diag/combat_diag_plot-2021-06-01.pdf`)
- normalized expression values (not included)
- spectra summary info (`./SpectraData/pca_details-2021-06-01.csv`)
- the rotation matrix and gene centers (`./SpectraData/spectra_rotation_matrix-gene_centers-2021-06-01.csv.gz`)
- spectra variables including the ComBat covariates (`./SpectraData/spectra_scores-batch_variables-2021-06-01.csv `)

### post-processing 

The spectra variables are centered by default and we divide each by its standard deviation to put each spectrum on the same scale. All subsequent modeling steps were performed using these centered and scaled spectra variables. 



## predictive modeling





#### bootstrap internal validation 





## descriptive modeling 





## added predictive value (APV)

https://www.fharrell.com/post/addvalue/



## follow-up samples 





