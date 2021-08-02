### CoMMpass preprocessing for tissue_spectra_process 20210601 
## using transcript expression 
## baseline samples only for derivation, but follow-up included as correct_only 
## protein coding genes, excluding ChrY, X and M, excluding Arora discordant genes
## ***this version DOES NOT include genes on the X chromosome*** 

### bring in expression data 
library(vroom)
## all transcript count data
tx_url <- "~/Documents/stats/survival/commpass_proj/commpass_spectra/MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz" 
tx_counts <- vroom(tx_url)
colnames(tx_counts)[1] <- "transcript_id"

## checks
# sum(is.na(tx_counts)) # result should be zero, can be slow

## peek at tx data
# dim(tx_counts) # 194715    909
# tx_counts[1:5,1:5]

### bring in annotation file 
gannot_url <- "~/Documents/stats/survival/commpass_proj/commpass_spectra/Homo_sapiens.GRCh37.74.gtf.flat.csv" 
gencode_annot <- vroom(gannot_url)

## peek at annotation data
# dim(gencode_annot) # 2244857      15
# head(gencode_annot)

library(readr)
library(dplyr)
library(tidyr)

### gene metadata 

## only need the gene_id and gene_name cols for this part 
gene_md <- gencode_annot %>% 
  distinct(gene_id, gene_name, .keep_all=T) %>% 
  select(-attributes)

## genes that need to "report" status for UAMS and SBUK
## AL035071.1 updated to RP5-1085F17.3 and SELENOI updated to EPT1 20210601 
## based on gene ids and matching  names in Homo_sapiens.GRCh37.74.gtf.flat.csv
UAMS <- c("EPT1","PDHA1","CHML","PFN1","YWHAZ","RAN","TAGLN2","ALDOA","CBX3","LGALS1",
          "ENO1","RUVBL1","CKS1B","CCT2","FABP5","TMPO","LARS2","RFC4","TRIP13","AURKA",
          "KIF20B","IFI16","KIF14","AIM2","LAS1L","ILF3","BIRC5","PSMD4","SLC19A1","WEE1",
          "ROBO1","AGO2","UBE2I","NADK","TBCB","DSG2","MTPAP","ASPM","TBRG4","NOL11","SNX5",
          "RAD18","CMSS1","CPSF3","FAM72B","CENPW","ZNF829","CHRM3","NXPE3","TCOF1","EXOSC4",
          "ARHGAP29","LINC01016","AHCYL1","GNG10","LTBP1","FUCA1","EVI5","PNPLA4",
          "TRIM33","CLCC1","CTBS","TMEM167B","ITPRIP","UBE2R2","TAF13","RP5-1085F17.3",
          "TRIM13","SLC48A1")

SBUK <- c("ADSS","BIRC5","CACYBP","CCT2","CCT3","CKS1B",
          "CTPS1","ENO1","GAPDH","KIAA0101","MSH6","NONO",
          "PRKDC","RAN","SF3B4","TFB2M","UBE2A")

## genes to remove from Arora discordant genes paper 
dqgenes <- 
  read_csv("~/Documents/camp_lab/rna_seq/Arora_DQ_2020/Union_Discordant_genes_TCGA_GTEx_Arora_TS2B.csv")
colnames(dqgenes)[1] <- "gene_name"

## CoMMpass note: the gene named "CKS1B" has 2 ENSG ids, info as of 20210524 from ensembl web...
# ENSG00000173207 actual protein coding gene (on chr1)
# ENSG00000268942 marked as pseudogene on chr5
## also note, the following UAMS and SBUK genes are on "patch" chromosomes and have duplicate
# gene_ids but only one is in expression counts
# SF3B4, FAM72B, EXOSC4 

ok_chroms <- c(as.character(1:22)) # X is excluded 

gene_md <- gene_md %>% 
  mutate(gene_status=case_when(gene_name %in% UAMS ~ "report",
                               gene_name %in% SBUK ~ "report", 
                               !chr %in% ok_chroms ~ "remove",
                               gene_name %in% dqgenes$gene_name ~ "remove",
                               gene_biotype=="protein_coding" ~ "allow",
                               TRUE ~ "remove")) 

# table(gene_md$gene_status, useNA="ifany") 
# allow remove report 
# 17324  46268     85 

## write out gene metadata file 
gene_md %>%
  select(gene_id, gene_name, gene_status) %>%
  filter(gene_status != "remove") %>%
  write_csv(., file="commpass_gene_metadata_20210601.csv")

### putting transcript length in count file

## if your tx_counts file doesn't contain gene_id, uncomment the following: 
tx_counts <- left_join(tx_counts, distinct(gencode_annot, transcript_id, gene_id), by="transcript_id") %>% 
  relocate(gene_id, transcript_id)

## get transcript length  
gencode_tcr_len <- gencode_annot %>% filter(feature=="exon") %>% 
  mutate(ex_size=(end-start)+1) %>% group_by(transcript_id) %>% 
  dplyr::summarise(length = sum(ex_size))

## join on transcript_id to add transcript lengths for genes remaining
tx_counts <- 
  inner_join(tx_counts, gencode_tcr_len, by="transcript_id") %>% 
  relocate(length, .after="transcript_id")

## there are 93 transcripts in the counts that have no transcript and no gene in the annotation
## with left_join the gene lengths and gene names are NAs 
## an inner_join replacing left_join removed the NAs and the few bogus trasncripts

## pulls rows with NAs
# tx_counts[which(rowSums(is.na(tx_counts))>0),] 

## write out transcript counts with transcript length and gene_id cols
vroom_write(tx_counts, "commpass_transcript_counts_20210601.gz",
            delim=",")

### sample metadata 

## per patient visit data for sample type (use baseline for derivation)
cmps_ppv <- 
  read_csv("~/Documents/stats/survival/commpass_proj/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv",
           guess_max=17321) %>% 
  select(SPECTRUM_SEQ, VJ_INTERVAL)
cmps_ppv <- cmps_ppv %>% drop_na(SPECTRUM_SEQ) 

## there are two duplicated SPECTRUM_SEQ values, but neither appears to have RNAseq
# cmps_ppv$SPECTRUM_SEQ[duplicated(cmps_ppv$SPECTRUM_SEQ)]
## mislabeled sample? again, doesn't have RNAseq 
# cmps_ppv %>% filter(SPECTRUM_SEQ %in% cmps_ppv$SPECTRUM_SEQ[duplicated(cmps_ppv$SPECTRUM_SEQ)])

## variables for combat  
# batch variable for combat
seq_qcsum <- 
  read_csv("~/Documents/stats/survival/commpass_proj/commpass_spectra/MMRF_CoMMpass_IA14_Seq_QC_Summary.csv",
           guess_max=5000) %>% 
  select(public_id=`Patients::KBase_Patient_ID`, visit_id=`Visits::Study Visit ID`, # reason=`Visits::Reason_For_Collection`,
         sample_id=`QC Link SampleName`, batch=Batch, 
         release_status=MMRF_Release_Status)
# fix sample_id col so compatible with other data (tx counts)
seq_qcsum <- seq_qcsum %>% 
  filter(release_status=="RNA-Yes") %>% 
  separate(sample_id, into=c("sid1","sid2","sid3","sid4",NA,NA,NA,NA)) %>% 
  unite(col=sample_id, sid1,sid2,sid3,sid4)
# age and gender covariates
cmps_perpat <- 
  read_csv("~/Documents/stats/survival/commpass_proj/commpass_spectra/MMRF_CoMMpass_IA14_PER_PATIENT.csv") %>% 
  select(public_id=PUBLIC_ID, D_PT_age, D_PT_gender)
# survival covariates 
cmps_sosurv <- 
  read_csv("~/Documents/stats/survival/commpass_proj/commpass_spectra/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv") %>% 
  select(public_id,ttcos,censos,ttcpfs,censpfs,ttctf1,censtf1)

# variable to distinguish baseline from follow-up samples 
# (simple 2 categories, VISITDY has quant info but many baseline values are negative)
cmps_ppv <- cmps_ppv %>% 
  mutate(smp_group=case_when(VJ_INTERVAL=="Baseline" ~ "Baseline", TRUE ~ "Follow-up"))

## combine variables to make sample metadata file
# get ids and batch from sequence QC 
# remove samples that are not from bone marrow (_BM)
# set status, only using baseline samples for spectra derivation 
sample_data <- left_join(seq_qcsum, cmps_ppv, by=c("visit_id"="SPECTRUM_SEQ")) %>% 
  select(-visit_id, -release_status) %>% 
  filter(grepl('BM', sample_id)) %>% 
  mutate(status=case_when(VJ_INTERVAL=="Baseline" ~ "use", 
                          TRUE ~ "correct_only")) %>% 
  drop_na() # removes any row with NA (there are 4 RNAseq samples with NA VJ_INTERVAL and NA survival covars)

# join age and gender (as factor)
sample_data <- left_join(sample_data, cmps_perpat, by="public_id")
sample_data$D_PT_gender <- factor(sample_data$D_PT_gender, levels=c(1,2), labels=c("male","female"))
# join survival covars
sample_data <- left_join(sample_data, cmps_sosurv, by="public_id")
# clean up and order cols
sample_data <- sample_data %>% 
  select(-public_id, -VJ_INTERVAL) %>% 
  relocate(status, .after="sample_id")

write_csv(sample_data, "commpass_sample_data_auto_followincl_20210601.csv")

