# CoMMpass MM spectra derivation and modeling 

## intro

Code used to generate spectra using the MMRF  CoMMpass multiple myeloma (MM) RNAseq data.  

more on CoMMpass: https://themmrf.org/finding-a-cure/our-work/the-mmrf-commpass-study/ 



## input data 

The current CoMMpass data is in the NCI GDC data portal but is controlled access through sbGaP. 

https://portal.gdc.cancer.gov/projects/MMRF-COMMPASS 

https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000748.v7.p4

We used versions IA14 of the RNAseq transcript counts, patient data, and sample QC data that were previously available from the MMRF. 

### additional reference data 

We also used a version of `Homo_sapiens.GRCh37.74.gtf` that was converted to flat text and saved as a csv file (deposited here as `Homo_sapiens.GRCh37.74.gtf.flat.csv`) 

We removed genes identified by Arora et al. (https://doi.org/10.1038/s41598-020-59516-z) as discordant using  the file `Union_Discordant_genes_TCGA_GTEx_Arora_TS2B.csv` 

## spectra derivation 





## predictive modeling





## descriptive modeling 
