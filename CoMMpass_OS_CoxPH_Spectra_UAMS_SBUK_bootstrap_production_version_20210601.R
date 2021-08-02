### CoMMpass OS CoxPH Spectra UAMS SBUK bootstrap production version 20210601
## using spectra derived from baseline samples only, other bone marrow samples corrected only 
## Nominal pval variable selection inside bootstraps  
## GMM2 threshold for CoMMpass spectra (UAMS and SBUK have their own methods)

library(vroom)
library(readr)
library(dplyr)
library(tidyr)
library(survival)
library(Hmisc) # for concordance calcs
library(mclust) # for GMM thresholds

## bring in spectra and survival outcomes 
commpass_os <- vroom("./CoMMpass_inclfollowup_final_SpectraData_20210601/spectra_scores-batch_variables-2021-06-01.csv")
colnames(commpass_os) <- gsub(x=colnames(commpass_os), pattern="PC",replacement="Sp")

# filter to only baseline samples and only select variables needed for this analysis
# divide each spectra var by SD 
cmps_os_wid <- commpass_os %>% 
  filter(smp_group=="Baseline") %>% 
  select(sample_id, ttcos, censos, Sp1:Sp39) %>% 
  mutate(across(Sp1:Sp39, ~ .x/sd(.x)))

## uncomment to remove unscaled version of spectra with follow-up samples to reduce memory use 
# rm(commpass_os)

## uncomment to write out SD scaled version of spectra with OS data for baseline samples only for modeling
# write_csv(cmps_os_wid, "commpass_spectra_baseline_os_20210601.csv")

## calculates actuarial hazards ratio from survdiff object 
act.haz <- function(data, time, event, class){
  data <- as.data.frame(data)
  time <- data[,time]
  event <- data[,event]
  class <- data[,class]
  temp_sdf <- survdiff(Surv(time, event) ~ class)
  HR <- (temp_sdf$obs[2]/temp_sdf$exp[2])/(temp_sdf$obs[1]/temp_sdf$exp[1])
  return(HR)
}

### set-up for UAMS and SBUK classifiers
# expression values for genes in these 2 classifiers, sample_id is same as spectra
exp_cbat <- vroom("./CoMMpass_inclfollowup_final_SpectraData_20210601/normalized_expression-2021-06-01.csv.gz")

gencode_annot <- vroom("~/Documents/stats/survival/commpass_proj/commpass_spectra/Homo_sapiens.GRCh37.74.gtf.flat.csv") %>% 
  select(gene_id, gene_name) %>% 
  distinct(gene_id, .keep_all=T)

## UAMS
## list of up regulated genes updated 20210407
## AL035071.1 updated to RP5-1085F17.3 and SELENOI updated to EPT1 20210601 
## based on gene ids and matching  names in Homo_sapiens.GRCh37.74.gtf.flat.csv
up <- c("EPT1","PDHA1","CHML","PFN1","YWHAZ","RAN","TAGLN2","ALDOA","CBX3","LGALS1",
        "ENO1","RUVBL1","CKS1B","CCT2","FABP5","TMPO","LARS2","RFC4","TRIP13","AURKA",
        "KIF20B","IFI16","KIF14","AIM2","LAS1L","ILF3","BIRC5","PSMD4","SLC19A1","WEE1",
        "ROBO1","AGO2","UBE2I","NADK","TBCB","DSG2","MTPAP","ASPM","TBRG4","NOL11","SNX5",
        "RAD18","CMSS1","CPSF3","FAM72B","CENPW","ZNF829","CHRM3","NXPE3","TCOF1","EXOSC4")
## list of down regulated genes updated 20210407
down <- c("ARHGAP29","LINC01016","AHCYL1","GNG10","LTBP1","FUCA1","EVI5","PNPLA4",
          "TRIM33","CLCC1","CTBS","TMEM167B","ITPRIP","UBE2R2","TAF13","RP5-1085F17.3",
          "TRIM13","SLC48A1")
## SBUK genes and beta coefficients 
bjh17 <- data.frame("genes" = c("ADSS","BIRC5","CACYBP","CCT2","CCT3","CKS1B",
                                "CTPS1","ENO1","GAPDH","KIAA0101","MSH6","NONO",
                                "PRKDC","RAN","SF3B4","TFB2M","UBE2A"),
                    "coeff" = c(0.67,0.39,0.56,1.15,0.61,0.59,0.59,0.84,0.49,
                                0.30,0.91,1.03,1.04,0.65,1.01,0.83,0.66))
## get gene_ids for genes in signatures above
sign_gene_ids <- gencode_annot %>% 
  filter(gene_name %in% up | gene_name %in% down | gene_name %in% bjh17$genes)

## the gene named "CKS1B" has 2 ENSG ids, info as of 20210524 from ensembl web...
## ENSG00000173207 actual protein coding gene (on chr1)
## ENSG00000268942 marked as pseudogene on chr5
## max normalized log2 expression of ENSG00000268942 is -4.075
# exp_cbat %>% filter(gene_id=="ENSG00000268942") %>% select(gene_id, adjlogcpkmed) %>% summary()
## removing this gene_id

sign_gene_ids <- sign_gene_ids %>% filter(gene_id!="ENSG00000268942")

## filter gene expression for genes in signatures above, only baseline samples 
## use gene names. genes are columns, individuals are rows
exp_cbat_w <- exp_cbat %>% 
  filter(gene_id %in% sign_gene_ids$gene_id & smp_group=="Baseline") %>% 
  left_join(., sign_gene_ids, by="gene_id") %>% 
  select(sample_id, gene_name, adjlogcpkmed) %>% 
  pivot_wider(names_from="gene_name", values_from="adjlogcpkmed") %>% 
  arrange(sample_id)

### calculating once on all data 
## UAMS
# select up regulated genes in data
anno_up <- exp_cbat_w %>% dplyr::select(intersect(up,colnames(exp_cbat_w))) 
# select down regulated genes in data
anno_dw <- exp_cbat_w %>% dplyr::select(intersect(down,colnames(exp_cbat_w))) 
## expression already in log2 scale
# compute means
uams_results <- exp_cbat_w[,"sample_id"] 
colnames(uams_results) <- "sample_id"
uams_results$up <- anno_up %>% rowMeans()
uams_results$dw <- anno_dw %>% rowMeans()
# compute proportion of mean up/down
uams_results$ua_score <- uams_results$up - uams_results$dw # expression already in log2 scale
## Find threshold by k-means clustering
# k-means cluster analysis
fit <- kmeans(uams_results$ua_score, 3)
# get cluster max
c.uams_all <- aggregate(uams_results$ua_score,by=list(fit$cluster),FUN=max)
# append cluster assignment
uams_results$fit_cluster <- fit$cluster
# classify
uams_results$ua_group = if_else(uams_results$ua_score > median(c.uams_all$x),"high","low")
uams_results$set <- "all_data"
# combine results
uams_results_sdata <- left_join(uams_results, 
                                dplyr::select(cmps_os_wid, sample_id, ttcos, censos), by="sample_id") 
uams_results_sdata$ua_group <- factor(uams_results_sdata$ua_group, levels = c("low","high"))
## actuarial HR
uams_act_hr <- act.haz(uams_results_sdata, "ttcos", "censos", "ua_group")

## SBUK
# get genes 
anno_17_all <- exp_cbat_w %>% dplyr::select(intersect(bjh17$genes,colnames(exp_cbat_w))) 
sbuk_results <- exp_cbat_w[,"sample_id"] 
# calculate the score with fixed coef from paper
sbuk_fixed_coeff <- bjh17 %>% 
  filter(genes %in% intersect(bjh17$genes,colnames(exp_cbat_w))) %>% 
  pull(coeff)
sbuk_results$sb_score <- as.numeric(data.matrix(anno_17_all) %*% sbuk_fixed_coeff)
# find cutoff 
high.cut_all <- quantile(sbuk_results$sb_score,.75) 
# classify
sbuk_results$sb_group <- if_else(sbuk_results$sb_score > high.cut_all,"high","low")
sbuk_results$set <- "all_data"
# combine results 
sbuk_results_sdata <- left_join(sbuk_results, 
                                dplyr::select(cmps_os_wid, sample_id, ttcos, censos), by="sample_id") 
sbuk_results_sdata$sb_group <- factor(sbuk_results_sdata$sb_group, levels = c("low","high"))
## actuarial HR
sbuk_act_hr <- act.haz(sbuk_results_sdata, "ttcos", "censos", "sb_group")

## write out results of these classifiers on all data
write_csv(bind_cols(uams_act_hr=uams_act_hr, sbuk_act_hr=sbuk_act_hr), 
          "commpass_os_alldata_uams_sbuk_HRs_moregenes_20210601.csv") # summary results for HRs from UAMS and SBUK

## bootstrapping 
set.seed(5197)
brs <- 1000 # 1000 bootstrap replications 
## all spectra dataset 
cmps_os <- cmps_os_wid %>% dplyr::select(-sample_id)
cmps_class <- cmps_os %>% dplyr::select(ttcos, censos)
## store output
# list of sample_ids for bootstrap sample
bootlist_results <- matrix(nrow=nrow(cmps_os), ncol=brs)
# C index all vars
cindex_b_boot <- vector("numeric", length=brs) # estimates of C on bootstrap samples using model output
cindex_b_orig_fh <- vector("numeric", length=brs) # estimates of C on original data
cindex_oob_fh <- vector("numeric", length=brs) # estimates of C on oob sample

# C index sigp (p<0.05) vars
cindex_b_boot_sigp <- vector("numeric", length=brs) # estimates of C on bootstrap samples using model output
cindex_b_orig_sigp <- vector("numeric", length=brs) # estimates of C on original data
cindex_oob_sigp <- vector("numeric", length=brs) # estimates of C on oob sample

# linear predictors from all vars CoxPH model (basis for classification)
cph_boot_linear_preds <- matrix(nrow=nrow(cmps_os), ncol=brs)
cph_orig_linear_preds <- matrix(nrow=nrow(cmps_os), ncol=brs)

# linear predictors from all sigp (p<0.05) vars CoxPH model (basis for classification)
cph_boot_linear_preds_sigp <- matrix(nrow=nrow(cmps_os), ncol=brs)
cph_orig_linear_preds_sigp <- matrix(nrow=nrow(cmps_os), ncol=brs)

# bootstrap outcome predictions for orig data (classes based on GMM threshold)
gmm2_orig_classes <- matrix(nrow=nrow(cmps_os), ncol=brs) 
gmm2_orig_classes_sigp <- matrix(nrow=nrow(cmps_os), ncol=brs)

# actuarial HRs 
bs_gmm2_act_hr <- vector("numeric", length=brs)
or_gmm2_act_hr <- vector("numeric", length=brs)
oob_gmm2_act_hr <- vector("numeric", length=brs)
bs_sigp_act_hr <- vector("numeric", length=brs)
or_sigp_act_hr <- vector("numeric", length=brs)
oob_sigp_act_hr <- vector("numeric", length=brs)

# thresholds
all_thresh_results <- vector("numeric", length=brs)
sigp_thresh_results <- vector("numeric", length=brs)

# UAMS and SBUK thresholds 
uams_thresholds <- vector(mode="numeric", length=brs+1)
uams_thresholds[1] <- median(c.uams_all$x)
sbuk_set_thresholds <- vector(mode="numeric", length=brs+1)
sbuk_set_thresholds[1] <- high.cut_all

# actuarial HRs for UAMS and SBUK 
uams_hr_boot <- vector("numeric", length=brs)
uams_hr_orig <- vector("numeric", length=brs)
uams_hr_oob <- vector("numeric", length=brs)
sbuk_hr_boot <- vector("numeric", length=brs)
sbuk_hr_orig <- vector("numeric", length=brs)
sbuk_hr_oob <- vector("numeric", length=brs)

# pvals
boot_pvals_alls <- data.frame(spectra=paste0("Sp",1:39))
boot_pvals_sigp <- data.frame(spectra=paste0("Sp",1:39))

for(b in 1:brs){
  ## using all vars
  boot_list <- sample(1:nrow(cmps_os), size = nrow(cmps_os), replace=T)
  bootlist_results[,b] <- cmps_os_wid$sample_id[boot_list]
  boot_set <- cmps_os[boot_list,]
  bs_class <- boot_set %>% dplyr::select(ttcos, censos)
  temp_mod <- coxph(Surv(ttcos, censos) ~ ., data=boot_set)
  boot_pvals_alls <- bind_cols(boot_pvals_alls, summary(temp_mod)$coefficients[,"Pr(>|z|)"], .name_repair="minimal")
  tmoslp <- temp_mod$linear.predictors
  cph_boot_linear_preds[,b] <- tmoslp
  bs_class$lp <- tmoslp
  cindex_b_boot[b] <- temp_mod$concordance["concordance"]
  # oob sample
  oob_set <- cmps_os[-boot_list,]
  oob_class <- oob_set %>% dplyr::select(ttcos, censos)
  oob_class$lp <- predict(temp_mod, oob_set %>% dplyr::select(contains("Sp")))
  cindex_oob_fh[b] <- 1-rcorr.cens(oob_class$lp, Surv(oob_class$ttcos, oob_class$censos))[1]
  # original data
  cmps_class$lp <- predict(temp_mod, cmps_os %>% dplyr::select(contains("Sp")))
  cph_orig_linear_preds[,b] <- cmps_class$lp
  cindex_b_orig_fh[b] <- 1-rcorr.cens(cmps_class$lp, Surv(cmps_class$ttcos, cmps_class$censos))[1]
  
  ## var reduced real bootstrap models 
  sigp_sps <- names(summary(temp_mod)$coefficients[,"Pr(>|z|)"][summary(temp_mod)$coefficients[,"Pr(>|z|)"]<0.05])
  sigp_mod <- coxph(Surv(ttcos, censos) ~ ., 
                    data=dplyr::select(boot_set, ttcos, censos, all_of(sigp_sps)))
  sigp_pval <- data.frame(pvals=summary(sigp_mod)$coefficients[,"Pr(>|z|)"]) %>% 
    tibble::rownames_to_column(var="spectra")
  boot_pvals_sigp <- left_join(boot_pvals_sigp, sigp_pval, by="spectra")
  sigp_tmoslp <- sigp_mod$linear.predictors
  cph_boot_linear_preds_sigp[,b] <- sigp_tmoslp
  bs_class$sigp_lp <- sigp_tmoslp
  cindex_b_boot_sigp[b] <- sigp_mod$concordance["concordance"]
  # oob sample
  oob_set <- cmps_os[-boot_list,]
  oob_class$sigp_lp <- predict(sigp_mod, oob_set %>% dplyr::select(all_of(sigp_sps)))
  cindex_oob_sigp[b] <- 1-rcorr.cens(oob_class$sigp_lp, Surv(oob_class$ttcos, oob_class$censos))[1]
  # original data
  cmps_class$sigp_lp <- predict(sigp_mod, cmps_os %>% dplyr::select(all_of(sigp_sps)))
  cph_orig_linear_preds_sigp[,b] <- cmps_class$sigp_lp
  cindex_b_orig_sigp[b] <- 1-rcorr.cens(cmps_class$sigp_lp, Surv(cmps_class$ttcos, cmps_class$censos))[1]
  
  ### GMM thresholds from mclust, only have actuarial and not CoxPH HR 
  ### G=2 threshold, equal variance, all vars
  ## on bootstrap data
  gmm_boot_fit2 <- densityMclust(bs_class$lp, G=2, modelNames="E", verbose=F)
  gmm_boot_fit2_xc_df <- data.frame(x=gmm_boot_fit2$data, 
                                    class=as.factor(gmm_boot_fit2$classification))
  all_thresh_results[b] <-  
    (gmm_boot_fit2_xc_df %>% arrange(x) %>% filter(class==2) %>% slice_head() %>% pull(x) + 
       gmm_boot_fit2_xc_df %>% arrange(x) %>% filter(class==1) %>% slice_tail() %>% pull(x))/2
  bs_class$gmm2_class <- ifelse(gmm_boot_fit2$classification==2, "high", "low")
  bs_class$gmm2_class <- factor(bs_class$gmm2_class, levels=c("low", "high"))
  ## actuarial HR
  bs_gmm2_act_hr[b] <- act.haz(bs_class, "ttcos", "censos", "gmm2_class")
  ## on original data
  cmps_class$cl_gmm2 <- max.col(predict(gmm_boot_fit2, newdata=cmps_class$lp, what="z"))
  cmps_class$gmm2_class <- ifelse(cmps_class$cl_gmm2==2, "high", "low")
  cmps_class$gmm2_class <- factor(cmps_class$gmm2_class, levels=c("low", "high"))
  gmm2_orig_classes[,b] <- cmps_class$gmm2_class
  ## actuarial HR
  or_gmm2_act_hr[b] <- act.haz(cmps_class, "ttcos", "censos", "gmm2_class")
  ## on oob sample 
  oob_class$cl_gmm2 <- max.col(predict(gmm_boot_fit2, newdata=oob_class$lp, what="z"))
  oob_class$gmm2_class <- ifelse(oob_class$cl_gmm2==2, "high", "low")
  oob_class$gmm2_class <- factor(oob_class$gmm2_class, levels=c("low", "high"))
  ## actuarial HR
  if(table(oob_class$gmm2_class)["low"]==0 | table(oob_class$gmm2_class)["high"]==0){
    oob_gmm2_act_hr[b] <- NA
  }else{
    oob_gmm2_act_hr[b] <- act.haz(oob_class, "ttcos", "censos", "gmm2_class")
  }
  
  ### G=2 threshold, equal variance, reduced vars sigp
  ## on bootstrap data
  gmm_boot_fit2_sigp <- densityMclust(bs_class$sigp_lp, G=2, modelNames="E", verbose=F)
  gmm_boot_fit2_xc_df_sigp <- data.frame(x=gmm_boot_fit2_sigp$data, 
                                         class=as.factor(gmm_boot_fit2_sigp$classification))
  sigp_thresh_results[b] <-  
    (gmm_boot_fit2_xc_df_sigp %>% arrange(x) %>% filter(class==2) %>% slice_head() %>% pull(x) + 
       gmm_boot_fit2_xc_df_sigp %>% arrange(x) %>% filter(class==1) %>% slice_tail() %>% pull(x))/2
  bs_class$sigp_class <- ifelse(gmm_boot_fit2_sigp$classification==2, "high", "low")
  bs_class$sigp_class <- factor(bs_class$sigp_class, levels=c("low", "high"))
  ## actuarial HR
  bs_sigp_act_hr[b] <- act.haz(bs_class, "ttcos", "censos", "sigp_class")
  ## on original data
  cmps_class$cl_sigp <- max.col(predict(gmm_boot_fit2_sigp, newdata=cmps_class$sigp_lp, what="z"))
  cmps_class$sigp_class <- ifelse(cmps_class$cl_sigp==2, "high", "low")
  cmps_class$sigp_class <- factor(cmps_class$sigp_class, levels=c("low", "high"))
  gmm2_orig_classes_sigp[,b] <- cmps_class$sigp_class
  ## actuarial HR
  or_sigp_act_hr[b] <- act.haz(cmps_class, "ttcos", "censos", "sigp_class")
  ## on oob sample 
  oob_class$cl_sigp <- max.col(predict(gmm_boot_fit2_sigp, newdata=oob_class$sigp_lp, what="z"))
  oob_class$sigp_class <- ifelse(oob_class$cl_sigp==2, "high", "low")
  oob_class$sigp_class <- factor(oob_class$sigp_class, levels=c("low", "high"))
  ## actuarial HR
  if(table(oob_class$sigp_class)["low"]==0 | table(oob_class$sigp_class)["high"]==0){
    oob_sigp_act_hr[b] <- NA
  }else{
    oob_sigp_act_hr[b] <- act.haz(oob_class, "ttcos", "censos", "sigp_class")
  }

  ## UAMS from bootstrap sample on bootstrap sample using k=3
  anno_up <- exp_cbat_w[boot_list,] %>% dplyr::select(intersect(up,colnames(exp_cbat_w))) 
  anno_dw <- exp_cbat_w[boot_list,] %>% dplyr::select(intersect(down,colnames(exp_cbat_w))) 
  uams <- exp_cbat_w[boot_list,"sample_id"]
  uams$up <- anno_up %>% rowMeans()
  uams$dw <- anno_dw %>% rowMeans()
  uams$ua_score <- uams$up - uams$dw 
  fit <- kmeans(uams$ua_score, 3)
  c.uams <- aggregate(uams$ua_score,by=list(fit$cluster),FUN=max)
  uams_thresholds[b+1] <- median(c.uams$x)
  uams$fit_cluster <- fit$cluster
  uams$ua_group <- if_else(uams$ua_score > median(c.uams$x),"high","low")
  uams$set <- paste0("boot",b)
  uams_sdata <- left_join(dplyr::select(uams, sample_id, ua_group), 
                          dplyr::select(cmps_os_wid, sample_id, ttcos, censos), by="sample_id") 
  uams_sdata$ua_group <- factor(uams_sdata$ua_group, levels = c("low","high"))
  ## UAMS from bootstrap sample on original data
  uams_results_sdata_temp <- uams_results_sdata
  uams_results_sdata_temp$ua_gr_orig <- if_else(uams_results_sdata_temp$ua_score > median(c.uams$x),"high","low")
  uams_results_sdata_temp$ua_gr_orig <- factor(uams_results_sdata_temp$ua_gr_orig, levels = c("low","high"))
  ## UAMS from bootstrap sample on oob sample 
  uams_results_sdata_oob <- uams_results_sdata[-boot_list,]
  uams_results_sdata_oob$ua_gr_oob <- if_else(uams_results_sdata_oob$ua_score > median(c.uams$x),"high","low")
  uams_results_sdata_oob$ua_gr_oob <- factor(uams_results_sdata_oob$ua_gr_oob, levels = c("low","high"))
  
  ## SBUK on bootstrap sample with fixed coeff from paper using top 0.25
  anno_17 <- exp_cbat_w[boot_list,] %>% dplyr::select(intersect(bjh17$genes,colnames(exp_cbat_w))) 
  sbuk <- exp_cbat_w[boot_list,"sample_id"]
  sbuk$sb_score <- as.numeric(data.matrix(anno_17) %*% sbuk_fixed_coeff) # score using betas from paper
  high.cut <- quantile(sbuk$sb_score,.75) 
  sbuk_set_thresholds[b+1] <- high.cut
  sbuk$sb_group <- if_else(sbuk$sb_score > high.cut,"high","low")
  sbuk$set <- paste0("boot",b)
  sbuk_sdata <- left_join(dplyr::select(sbuk, sample_id, sb_group), 
                          dplyr::select(cmps_os_wid, sample_id, ttcos, censos), by="sample_id") 
  sbuk_sdata$sb_group <- factor(sbuk_sdata$sb_group, levels = c("low","high"))
  ## SBUK fixed beta from bootstrap sample on original data
  sbuk_results_sdata_temp <- sbuk_results_sdata
  sbuk_results_sdata_temp$sb_gr_orig <- if_else(sbuk_results_sdata_temp$sb_score > high.cut,"high","low")
  sbuk_results_sdata_temp$sb_gr_orig <- factor(sbuk_results_sdata_temp$sb_gr_orig, levels = c("low","high"))
  ## SBUK fixed beta from bootstrap sample on oob sample 
  sbuk_results_sdata_oob <- sbuk_results_sdata[-boot_list,]
  sbuk_results_sdata_oob$sb_gr_oob <- if_else(sbuk_results_sdata_oob$sb_score > high.cut,"high","low")
  sbuk_results_sdata_oob$sb_gr_oob <- factor(sbuk_results_sdata_oob$sb_gr_oob, levels = c("low","high"))

  ## UAMS
  # HRs on bootstrap data
  uams_hr_boot[b] <- act.haz(uams_sdata, "ttcos", "censos", "ua_group")
  # HRs on original data 
  uams_hr_orig[b] <- act.haz(uams_results_sdata_temp, "ttcos", "censos", "ua_gr_orig")
  # HRs on oob data
  uams_hr_oob[b] <- act.haz(uams_results_sdata_oob, "ttcos", "censos", "ua_gr_oob")

  ## SBUK with fixed betas
  # HRs on bootstrap data 
  sbuk_hr_boot[b] <- act.haz(sbuk_sdata, "ttcos", "censos", "sb_group")
  # HRs on original data 
  sbuk_hr_orig[b] <- act.haz(sbuk_results_sdata_temp, "ttcos", "censos", "sb_gr_orig")
  # HRs on oob data
  sbuk_hr_oob[b] <- act.haz(sbuk_results_sdata_oob, "ttcos", "censos", "sb_gr_oob")
}

bootlist_results <- data.frame(bootlist_results)
colnames(bootlist_results) <- paste0("boot", 1:brs)
spectra_actuarial_HR_res <- data.frame(bs_alls=bs_gmm2_act_hr, or_alls=or_gmm2_act_hr,
                                       oob_alls=oob_gmm2_act_hr, 
                                       bs_sigp=bs_sigp_act_hr, or_sigp=or_sigp_act_hr,
                                       oob_sigp=oob_sigp_act_hr) 
cindex_results <- data.frame(cindex_b_boot_alls=cindex_b_boot, cindex_b_orig_alls=cindex_b_orig_fh, 
                             cindex_oob_alls=cindex_oob_fh, 
                             cindex_b_boot_sigp, cindex_b_orig_sigp, 
                             cindex_oob_sigp)
cph_boot_linear_preds_alls <- data.frame(cph_boot_linear_preds)
colnames(cph_boot_linear_preds_alls) <- paste0("boot", 1:brs)
cph_orig_linear_preds_alls <- data.frame(sample_id=cmps_os_wid$sample_id, cph_orig_linear_preds)
colnames(cph_orig_linear_preds_alls) <- c("sample_id", paste0("boot", 1:brs))

cph_boot_linear_preds_sigp <- data.frame(cph_boot_linear_preds_sigp)
colnames(cph_boot_linear_preds_sigp) <- paste0("boot", 1:brs)
cph_orig_linear_preds_sigp <- data.frame(sample_id=cmps_os_wid$sample_id, cph_orig_linear_preds_sigp)
colnames(cph_orig_linear_preds_sigp) <- c("sample_id", paste0("boot", 1:brs))

colnames(boot_pvals_alls) <- c("spectra", paste0("boot", 1:brs))
colnames(boot_pvals_sigp) <- c("spectra", paste0("boot", 1:brs))

gmm2_orig_classes <- data.frame(sample_id=cmps_os_wid$sample_id, gmm2_orig_classes)
colnames(gmm2_orig_classes) <- c("sample_id", paste0("boot", 1:brs))
gmm2_orig_classes_sigp <- data.frame(sample_id=cmps_os_wid$sample_id, gmm2_orig_classes_sigp)
colnames(gmm2_orig_classes_sigp) <- c("sample_id", paste0("boot", 1:brs))

all_thresh_results <- data.frame(alls=all_thresh_results, sigp=sigp_thresh_results)

uams_sbuk_def_thresh <- data.frame(uams_km=uams_thresholds, sbuk_75=sbuk_set_thresholds)

uams_sbuk_act_HR <- data.frame(uams_hr_boot, uams_hr_orig, uams_hr_oob, 
                               sbuk_hr_boot, sbuk_hr_orig, sbuk_hr_oob) 

write_csv(bootlist_results, "commpass_os_bootlist_results_20210601.csv")
write_csv(cindex_results, "commpass_os_cindex_results_20210601.csv")

write_csv(cph_boot_linear_preds_alls, "commpass_os_cph_boot_linear_preds_alls_20210601.csv")
write_csv(cph_orig_linear_preds_alls, "commpass_os_cph_orig_linear_preds_alls_20210601.csv")
write_csv(cph_boot_linear_preds_sigp, "commpass_os_cph_boot_linear_preds_sigp_20210601.csv")
write_csv(cph_orig_linear_preds_sigp, "commpass_os_cph_orig_linear_preds_sigp_20210601.csv")

write_csv(boot_pvals_alls, "commpass_os_boot_pvals_alls_20210601.csv")
write_csv(boot_pvals_sigp, "commpass_os_boot_pvals_sigp_20210601.csv")

write_csv(all_thresh_results, "commpass_os_all_thresh_results_20210601.csv")
write_csv(uams_sbuk_def_thresh, "commpass_uams_sbuk_thresh_results_20210601.csv")

write_csv(gmm2_orig_classes, "commpass_os_gmm2_alls_orig_classes_20210601.csv")
write_csv(gmm2_orig_classes_sigp, "commpass_os_gmm2_sigp_orig_classes_20210601.csv")

write_csv(spectra_actuarial_HR_res, "commpass_os_spectra_actuarial_HR_res_20210601.csv")

write_csv(uams_sbuk_act_HR, "commpass_os_uams_sbuk_act_HR_res_mg_20210601.csv") 


