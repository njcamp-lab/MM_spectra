cat('\nLoading required libraries\n')
suppressPackageStartupMessages({
library(ggplot2)
library(data.table)
library(gridExtra)
library(foreach)
library(RColorBrewer)
library(sva)
})

# Constants --
TRANSCRIPT_MODE<-FALSE
BATCH_CORRECT<-FALSE
LOW_GENE_COUNT_THRESHOLD <- 100
LOW_GENE_PERCENTILE <- 0.05
SAMPLE_LOW_GENE_THRESHOLD <- 0.10 # what percentage of total good genes in a sample have to have LOW_GENE_COUNT_TRESHOLD counts
OUTPUT_DIR <- 'SpectraData'
DIAG_DIR <- 'SpectraData/Diag'

# Functions---- 
# Find inflection point (elbow) in scree plot of PC variances
elbow_finder<-function(data) {
  elbow.dt<-data[order(-value)][,idx:=.I-1]
  elbow.dt[,selected:=as.factor('Not used')]
  
  slope<-(min(elbow.dt$value)-max(elbow.dt$value))/(nrow(elbow.dt)-1)
  perpslope<-(-1/slope)
  intercept<-max(elbow.dt$value) # Is this accurate?
  elbow.dt[,perpcept:=value - perpslope*idx]
  
  elbow.dt[,y:=(perpcept*slope - intercept*perpslope)/(slope-perpslope)]
  elbow.dt[,x:=(intercept-perpcept)/(perpslope - slope)]
  elbow.dt[,dist:=sqrt((value-y)^2 + (idx-x)^2)]
  
  maxidx<-which.max(elbow.dt$dist)
  elbow.dt[idx<maxidx]$selected<-as.factor('Selected')
  elbow.dt[,propvar:=value/sum(value)]
  elbow.dt[,cumulpropvar:=cumsum(value)/sum(value)]
  elbow.dt[]
}

# Visualize elbow plot
# input data is output of elbow_finder
plot_elbow<-function(data) {
  slope<-(min(data$value)-max(data$value))/(nrow(data)-1)
  intercept<-max(data$value)
  
  g1<-ggplot(data)+
    geom_point(aes(y=value,x=idx,color=selected))+
    geom_abline(slope=slope,intercept=intercept,color='red')+
    theme_classic()+
    scale_color_manual(values=c('black','red'))+
    scale_shape_manual(values=c(16,15))+
    # geom_abline(intercept=intercept,slope=slope,linetype="dashed")+
    theme(legend.position='None',axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    ylab('Variance explained')+
    xlab('Spectra (PC)')+
    #scale_x_continuous(limits=c(-5,pc_count))+
    #geom_label_repel(data=pcvars.df[maxidx,],aes(x=idx,y=var,label=maxidx))+
    theme(axis.line.y=element_blank())+
    labs(title='Scree plot of variance explained per PC')
  g1
}


# for each of covariates that exists as a column in pca.summary, and for each column labeled "PC##" in pca.summary
# perform simple linear model (pc) and record F-test p-value
#  Inputs:
#   pca.summary: summary of PC variances and whether they are selected or not
#   pca.data: PC scores and other data for samples
#   covariates: names of covariates to test
pc.covariate.summary<-function(pca.data, pca.summary, covariates) {
  pcs<-colnames(pca.data)[colnames(pca.data)%like%'PC']
  vars<-colnames(pca.data)[colnames(pca.data)%in%covariates]

  # Sort PCs by number
  ftest.dt<-foreach(pc=pcs,.combine=rbind) %:%
    foreach(var=vars,.combine=rbind) %do% {
      lm<-lm(as.formula(paste(pc,'~',var)),pca.data)
      f<-summary(lm)$fstatistic
      data.table(component=pc,pc.num=as.numeric(gsub('PC','',pc)),variable=var,p=1-pf(f[1],f[2],f[3]))
    }
  # Reorder PC factor levels
  factor.levels<-data.table(level=levels(as.factor(ftest.dt$component)),levelnum=as.numeric(gsub('PC','',levels(as.factor(ftest.dt$component)))))
  ftest.dt$component<-factor(ftest.dt$component,levels=factor.levels[order(levelnum)]$level)
  
  merge(ftest.dt,pca.summary[,c('pc','propvar')],by.x='component',by.y='pc')
}

pc.summary<-function(data,batch.dt) {
    temp1.dt<-dcast(data,sample_id~gene_id,value.var='adjlogcpkmed')
    pca<-prcomp(temp1.dt[,.SD,.SDcols=!c('sample_id')],center=TRUE,scale=FALSE,retx=TRUE)
    pcvar<-data.table(pc=colnames(pca$x),value=pca$sdev^2)
    pca.summary.dt<-elbow_finder(pcvar)
    selected.pc.count<-nrow(pca.summary.dt[selected=='Selected'])
    pc.scores<-cbind(temp1.dt[,1],pca$x[,pca.summary.dt[selected=='Selected']$pc])
    pca.dt<-merge(batch.dt,pc.scores,by='sample_id')
    list(pca=pca,pca.summary=pca.summary.dt,data=pca.dt[])
  }
  

# Plot heatmap for phase-corrected PCs
# input is molten data.table output from pc.covariate.summary, with columns
#   component: PC as factor
#   pc.num: PC as number
#   variable: covariate we are testing
#   p: F-test p-value for linear model pc~variable
plot.pc.correction.heatmap<-function(molten.dt,midpoint=-10,subtitle='',levelnames='') {
  molten.dt[p<10^-20]$p<-10^-20 # Cap low p-values
  
  if(length(levelnames)>0) {
    molten.dt$variable<-factor(molten.dt$variable,levels=levelnames)
  } else {
    molten.dt$variable<-factor(molten.dt$variable)
  }
  pc.labels<-dcast(molten.dt,component+pc.num+propvar~1,fun.aggregate=max,value.var='propvar')[order(pc.num)]
  ggplot(molten.dt)+
    geom_tile(aes(x=pc.num,y=variable,fill=log10(p)),size=0.25,color='gray')+
    scale_fill_gradient2(low='red',mid='yellow',high='white',midpoint=midpoint,name="log(p)",limits=c(-20,0))+
    scale_x_continuous(name='Diagnostic PC',breaks=seq(1,nrow(pc.labels),by=1),labels=pc.labels$pc.num,
                       sec.axis=dup_axis(name='Variance',labels=sprintf("%0.3f",round(pc.labels$propvar,digits=3))))+
    theme_minimal()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())+
    ggtitle('Diagnostic PC association with covariates',subtitle=subtitle)
}


args = commandArgs(trailingOnly=TRUE)

# 1st Argument one is a file of expression counts, either transcript X sample or gene X sample; 2nd column should be gene/transcript length
# 2nd argument is the filename of the list of gene names to include; this is the 'good genes' list:
### It contains 3 columns: gene_id, gene name, and whether to allow it in, force it in, or exclude it 
### from the size-factor calculation and spectra derivation, but still calculate its normalized expression
# 3rd argument is the filename of the batch variables list; 1st column should be 'sample_id'; second column
### is either: 'use' or 'correct_only'; if 'use', this sample is used in the derivation of spectra
### if 'correct_only', then it is included in the normalization etc. but not derivation of spectra,
### third column is the batch variable; subsequent variables are protected in batch correction

# Main section-----
setwd(".")
if(length(args)<3) {
  stop('This script requires 3 arguments, in order: \n\t1. Expression counts file name\n\t2. Gene list filename\n\t3. Batch variable filename')
}


if(!dir.exists(OUTPUT_DIR)) {
  cat('\nCreating output directory ',OUTPUT_DIR)
  dir.create(OUTPUT_DIR)
} else {
  cat('\nDirectory',OUTPUT_DIR,'exists, data will be written there')
}


if(!dir.exists(DIAG_DIR)) {
  cat('\nCreating diagnostic directory ',DIAG_DIR)
  dir.create(DIAG_DIR)
} else {
  cat('\nDirectory',DIAG_DIR,'exists, diagnostics will be written there')
}

# Comment this out later 
expr.filename<-"../BLCA_BRCA_PRAD_COAD/combined_expression.csv.gz"
genes.filename<-'../tcga_genes.csv'
batch.filename<-'batch_brca_vars.csv'

expr.filename<-args[1]
genes.filename<-args[2]
batch.filename<-args[3]


# Load files, report on contents
expr.dt<-fread(expr.filename)
cat('\n\nLoaded expression datafile:',expr.filename)
cat('\n\t',nrow(expr.dt),'rows (transcripts/genes)')
cat('\n\t',ncol(expr.dt),'columns (samples)')

# Test here for transcript flavor v. gene flavor
TRANSCRIPT_MODE<-grepl('transcr',paste(colnames(expr.dt),collapse=' '))

if(TRANSCRIPT_MODE) {
  cat('\nTranscript column detected in expression file, running pipeline in transcript mode')
} else {
  cat('\nRunning pipeline in gene mode')
}

# gene_status options
# allow: genes are considered until they fail QC
# force: genes are evaluated by gene QC; if they pass, treat as normal; if they fail, exclude from sample QC but keep in rest of steps
# report: genes are evaluated by gene QC; if they pass, treat as normal; if they fail, exclude from sample QC, size-factor calc, spectra derivation

genes.dt<-fread(genes.filename)
cat('\n\nLoaded gene inclusion file:',genes.filename)
cat('\n\t',nrow(genes.dt),'genes loaded')
cat('\n\t',nrow(genes.dt[gene_status=='allow']),'genes allowed in size-factor calc and spectra derivation (prior to QC)')
force.genes<-genes.dt[gene_status=='force']
cat('\n\t',nrow(force.genes),'gene(s) forced into normalization (ignore QC) and spectra derivation')
report.genes<-genes.dt[gene_status=='report']
cat('\n\t',nrow(report.genes),'gene(s) allowed to bypass QC for reporting purposes only (corrected & normalized)')



# use: normal sample, use in every step unless it fails sample QC
# correct_only: unless it fails sample QC, use only in ComBat step; exclude from spectra derivation
batch.dt<-fread(batch.filename)
cat('\n\nLoaded batch variables file:',batch.filename)
cat('\n\t',nrow(batch.dt),'samples loaded')
cat('\n\t',nrow(batch.dt[status=='use']),'samples to be used in spectra derivation')
cat('\n\t',nrow(batch.dt[status=='correct_only']),'samples to be used only in batch correction')
cat('\n\t',ncol(batch.dt),'variables found (including sample_id & status)')

if(length(colnames(batch.dt))>2) {
  BATCH_CORRECT<-TRUE
  batch_count<-length(table(batch.dt[,3]))
  cat('\n\t',batch_count,'batch(es) found.')
  cat('\n\tSmallest batch size:',min(table(batch.dt[,3])))
}

# Melt expression data table
if(TRANSCRIPT_MODE) {
  temp.molten<-melt(expr.dt,id.vars=c('gene_id','transcript_id','length'),variable.name='sample_id',value.name='expected_count')
  molten.counts<-temp.molten[sample_id%in%batch.dt$sample_id]
  rm(temp.molten)
} else {
  temp.molten<-melt(expr.dt,id.vars=c('gene_id','length'),variable.name='sample_id',value.name='expected_count')
  molten.counts<-temp.molten[sample_id%in%batch.dt$sample_id]
  rm(temp.molten)
}
# Strip out any low-count genes
# Do this based on all included genes
low.genes<-molten.counts[sample_id%in%batch.dt[status=='use']$sample_id & gene_id%in%genes.dt$gene_id,
                         list(gene_count=sum(expected_count)),by=c('gene_id','sample_id')][gene_count<LOW_GENE_COUNT_THRESHOLD,.N,by='gene_id'][N/nrow(batch.dt[status=='use'])>=LOW_GENE_PERCENTILE]
# low genes includes forced and report only genes
cat('\n\nFiltered out',nrow(low.genes),'genes where more than',LOW_GENE_PERCENTILE,'of samples had fewer than',LOW_GENE_COUNT_THRESHOLD,'expected counts')

good.force.genes<-genes.dt[gene_status=='force' & !gene_id%in%low.genes$gene_id]
if(nrow(good.force.genes)>1) {
  cat('\n\n',nrow(good.force.genes),'genes were marked as "force" but pass QC. Updating them to "allow"')
  genes.dt[gene_id%in%good.force.genes$gene_id,gene_status:='allow']
} else if(nrow(good.force.genes)==1) {
  cat('\n\n',good.force.genes[1,]$gene_id,'was marked as "force", but passes QC. Updating it to "allow"')
  genes.dt[gene_id%in%good.force.genes$gene_id,gene_status:='allow']
}

good.report.genes<-genes.dt[gene_status=='report' & !gene_id%in%low.genes$gene_id]
if(nrow(good.report.genes)>1) {
  cat('\n\n',nrow(good.report.genes),'genes were marked as "report" but pass QC. Updating them to "allow"')
  genes.dt[gene_id%in%good.report.genes$gene_id,gene_status:='allow']
} else if (nrow(good.report.genes==1)) {
  cat('\n\n',good.report.genes[1,]$gene_id,'was marked as "report", but passes QC. Updating it to "allow"')
  genes.dt[gene_id%in%good.report.genes$gene_id,gene_status:='allow']
}

# These are the genes we will use to make sample QC decisions
# and derive size factors and eventually spectra
good.genes<-genes.dt[gene_id%in%molten.counts$gene_id & 
                       (gene_id%in%genes.dt[gene_status%in%c('allow')]$gene_id & !gene_id%in%low.genes$gene_id) | gene_id%in%force.genes$gene_id]
cat('\n\t',nrow(good.genes),'gene(s) remain')

# Look for poor samples; mark them as 'subpar' or something
cat('\n\tLooking for samples where more than 10% of the high-quality genes have fewer than',LOW_GENE_COUNT_THRESHOLD,'counts')
low.quality.samples<-molten.counts[gene_id%in%good.genes$gene_id & sample_id%in%batch.dt$sample_id,list(gene_count=sum(expected_count)),by=c('gene_id','sample_id')][gene_count<LOW_GENE_COUNT_THRESHOLD,.N,by=c('sample_id')][N>SAMPLE_LOW_GENE_THRESHOLD*nrow(good.genes)]
lowqual.use.count<-nrow(batch.dt[status=='use' & sample_id%in%low.quality.samples$sample_id])
lowqual.correct.count<-nrow(batch.dt[status=='correct_only' & sample_id%in%low.quality.samples$sample_id])

if(lowqual.correct.count>0) {
  cat('\n\t',lowqual.correct.count,'sample(s) marked as "correct_only" had low quality and will be removed.')
}

if(lowqual.use.count>0) {
  cat('\n\t',lowqual.use.count,'sample(s) marked as "use" had low quality and will be removed.')
}

batch.dt[sample_id%in%low.quality.samples$sample_id,status:='lowqual']
if(nrow(low.quality.samples)>0) {
  cat('\n\tIn total, removed',nrow(low.quality.samples),'low quality samples.')
}


# Process transcripts to per-gene counts
cat('\n\nProcessing transcript data to generate log2-transformed, median-normalized counts per kilobase')
temp.dt<-molten.counts[gene_id%in%good.genes$gene_id | gene_id%in%report.genes$gene_id]

cat('\nConverting gene/transcript length to kilobases')
temp.dt[,kb_length:=length/1000]

cat('\nMerging in gene status')
temp2.dt<-merge(temp.dt,genes.dt[gene_id%in%good.genes$gene_id | gene_id%in%report.genes$gene_id,c('gene_id','gene_status')],by='gene_id')

cat('\nGenerating counts/kb, summing across transcripts per sample/gene')
temp3.dt<-temp2.dt[,list(cpk=sum((expected_count+1/.N)/kb_length)),by=c('sample_id','gene_id','gene_status')]
cat('\nCalculating per-sample size factor as median counts/kb value')
size.factors<-temp3.dt[gene_id%in%good.genes$gene_id,list(size_factor=median(cpk)),by='sample_id']
final.dt<-merge(size.factors,temp3.dt,by='sample_id')
  
cat('\nScaling by size factor')
final.dt[,cpkmed:=cpk/size_factor]
cat('\nLog2-transforming scaled counts/kb')
final.dt[,logcpkmed:=log2(cpkmed)]

# Outlier truncation
cat('\nLooking for expression outliers (normalized value more than 5 SD from mean per gene)')
final.dt[,mean:=mean(logcpkmed),by='gene_id']
final.dt[,sd:=sd(logcpkmed),by='gene_id']
final.dt[,adjlogcpkmed:=logcpkmed]

outlier.count<-nrow(final.dt[(logcpkmed-mean)/sd>=5 | (logcpkmed-mean)/sd<= -5])
gene.outlier.count<-length(unique(final.dt[(logcpkmed-mean)/sd>=5 | (logcpkmed-mean)/sd<= -5]$gene_id))
sample.outlier.count<-length(unique(final.dt[(logcpkmed-mean)/sd>=5 | (logcpkmed-mean)/sd<= -5]$sample_id))
multi.outlier.genes<-final.dt[(logcpkmed-mean)/sd>=5 | (logcpkmed-mean)/sd<= -5,.N,by='gene_id'][N>1]  

cat('\n\tFound',outlier.count,'total outlier(s)')
cat('\n\t',gene.outlier.count,'gene(s) had at least one outlier')
cat('\n\t',nrow(multi.outlier.genes),'gene(s) had more than one outlier')
if(nrow(multi.outlier.genes>0)) {
  cat('\n\tMax outliers in a gene:',max(multi.outlier.genes$N))
  cat('\n\t\tSeen in',nrow(multi.outlier.genes[N==max(multi.outlier.genes$N)]),'gene(s)')
}
cat('\n\tTruncating outliers to mean +/- 5SD')
final.dt[(logcpkmed-mean)/sd>=5,adjlogcpkmed:=mean+5*sd]
final.dt[(logcpkmed-mean)/sd<= -5,adjlogcpkmed:=mean-5*sd]

# ComBat Batch Correction (on ALL samples)
  # get batch from argument, merge by sample_id to be sure of order
  # get covariates fom argument, merge by sample_id to be sure of order
  # dcast expr.molten to genes x sample matrix
cat('\nBatch correction using ComBat')

var.dt<-merge(dcast(final.dt[sample_id%in%batch.dt[status!='lowqual']$sample_id],sample_id+size_factor~1,value.var='adjlogcpkmed',fun.aggregate=length)[order(sample_id)][,1:3],batch.dt,by='sample_id')
if(BATCH_CORRECT) {
  combat.dt<-dcast(final.dt[sample_id%in%batch.dt[status!='lowqual']$sample_id],gene_id~sample_id,value.var='adjlogcpkmed')
  setcolorder(combat.dt,c('gene_id',batch.dt[status!='lowqual'][order(sample_id)]$sample_id))
  combat.matrix<-data.matrix(data.frame(combat.dt,row.names='gene_id',check.names=FALSE))
  
  
  batch<-as.factor(batch.dt[status!='lowqual'][order(sample_id)][[3]])
  batch.var.name<-colnames(batch.dt)[3]
  cat('\n\tUsing',batch.var.name,'as batch variable')

  # Check for NAs here?
  if(length(colnames(batch.dt))>3) {
    cat('\n\t','Protecting',length(colnames(batch.dt))-3,'covariate(s)')
    cat('\nRunning ComBat\n\n')
    cat('\n-------------Output from ComBat starts----------------------\n')
    cov.model<-model.matrix(as.formula(paste0('~',paste(colnames(batch.dt)[4:length(colnames(batch.dt))],collapse='+'))),batch.dt[status!='lowqual'][order(sample_id)])
    #cbat=ComBat(data.matrix(combat.dt[,-1]),batch,cov.model)
    cbat=ComBat(combat.matrix,batch,cov.model)
    cat('\n-------------Output from ComBat ends----------------------\n')
  } else {
    cat('\nRunning ComBat\n\n')
    cat('\n-------------Output from ComBat starts----------------------\n')
    cbat=ComBat(combat.matrix,batch)
    cat('\n-------------Output from ComBat ends----------------------')
  }
  cbat.dt=data.table(t(cbat),keep.rownames=TRUE)
  colnames(cbat.dt)[1]<-'sample_id'
  cat('\n\nCombat finished')
  #colnames(cbat.dt)<-as.character(combat.dt$gene_id)
  cat('\n\tMerging combat-corrected expression data with sample variables')
  # Need to get gene IDS back here
  #temp.molten<-melt(cbind(var.dt,cbat.dt),
  #                       id.vars=colnames(var.dt),
  #                       value.name='adjlogcpkmed',
  #                       variable.name='gene_id')
  temp.molten<-melt(merge(var.dt,cbat.dt,by='sample_id'),
                    id.vars=colnames(var.dt),
                    value.name='adjlogcpkmed',
                    variable.name='gene_id')
  final.molten<-merge(genes.dt[,c('gene_id','gene_status')],temp.molten,by='gene_id')
  uncorrected.molten<-merge(batch.dt[,c('sample_id','status',batch.var.name),with=FALSE],
                            final.dt[,c('sample_id','size_factor','gene_id','gene_status','adjlogcpkmed')],
                            by='sample_id')
  cor.compare.dt<-merge(final.molten[,c('sample_id','gene_id','adjlogcpkmed')],uncorrected.molten[,c('sample_id','gene_id','adjlogcpkmed',batch.var.name),with=FALSE],by=c('sample_id','gene_id'))
  cat('\n\tCalculating correlations between corrected and uncorrected data by gene')
  gene.correlations<-cor.compare.dt[,list(pearson_cor=cor(adjlogcpkmed.x,adjlogcpkmed.y)),by='gene_id']
  low_cor_genes<-nrow(gene.correlations[pearson_cor<0.9])
  cat('\n\t',low_cor_genes,'gene(s) found with Pearson correlation coefficient < 0.9')  
  if(mean(gene.correlations$pearson_cor)>0.999) {
    cat('\n\t***Gene correlations are unusually high, did ComBat do anything?***')
  }
  
  cat('\n\tCalculating correlations between corrected and uncorrected data by sample')
  sample.correlations<-cor.compare.dt[,list(pearson_cor=cor(adjlogcpkmed.x,adjlogcpkmed.y)),by='sample_id']
  low_cor_samples<-nrow(sample.correlations[pearson_cor<0.9])
  cat('\n\t',low_cor_samples,'sample(s) found with Pearson correlation coefficient < 0.9')  
  # Alert user here if these are too good or too bad?
  if(mean(sample.correlations$pearson_cor)>0.999) {
    cat('\n\t***Sample correlations are unusually high, did ComBat do anything?***')
  }
  
  gene_cor_filename=paste0(DIAG_DIR,'/correction_gene_correlation-',Sys.Date(),'.csv')
  fwrite(gene.correlations,gene_cor_filename)
  cat('\n\tCorrected/uncorrected gene correlations written to:',gene_cor_filename)
  
  sample_cor_filename=paste0(DIAG_DIR,'/correction_sample_correlation-',Sys.Date(),'.csv')
  fwrite(sample.correlations,sample_cor_filename)
  cat('\n\tCorrected/uncorrected samples correlations written to:',sample_cor_filename)
  
  # Diagnostic PCA & association tests
  cat('\n\tGenerating PCA for ComBat diagnostics')
  combat.pca<-pc.summary(final.molten[gene_status!='report'],batch.dt)
  nocombat.pca<-pc.summary(uncorrected.molten[gene_status!='report'],batch.dt)
  cat('\n\tComparing spectra (PCs) with sample variables for corrected data')
  combat.pca.cov.dt<-pc.covariate.summary(combat.pca$data,combat.pca$pca.summary,colnames(batch.dt)[-c(1:2)])
  combat.assoc.plot<-plot.pc.correction.heatmap(combat.pca.cov.dt,midpoint=-10,subtitle='Batch-corrected',levelnames=colnames(batch.dt)[-c(1:2)])
  
  cat('\n\tComparing spectra (PCs) with sample variables for un-corrected data')
  nocombat.pca.cov.dt<-pc.covariate.summary(nocombat.pca$data,nocombat.pca$pca.summary,colnames(batch.dt)[-c(1:2)])
  nocombat.assoc.plot<-plot.pc.correction.heatmap(nocombat.pca.cov.dt,midpoint=-10,subtitle='No Batch Correction',levelnames=colnames(batch.dt)[-c(1:2)])
  combat_diag_filename=paste0(DIAG_DIR,'/combat_diag_plot-',Sys.Date(),'.pdf')
  ggsave(combat_diag_filename,arrangeGrob(nocombat.assoc.plot,combat.assoc.plot),height=12,width=12)
  cat('\n\tCombat diagnostic plot written to:',combat_diag_filename)
} else {
  cat('\n\tSkipping batch correction.')
  final.molten<-merge(batch.dt[,c('sample_id','status')],final.dt[,c('sample_id','size_factor','gene_id','gene_status','adjlogcpkmed')],by='sample_id')
}

# Derive PCA on all genes not marked 'report' and all samples marked 'use'
spectra.pca<-pc.summary(final.molten[!gene_status%in%c('report') & status=='use'],batch.dt)
cat('\n\tFound',nrow(spectra.pca$pca.summary),'principal components')
cat('\n\tUsing elbow method to select PCs')
selected.pc.count<-nrow(spectra.pca$pca.summary[selected=='Selected'])
cat('\n\t\tSelected',selected.pc.count,
    'principal components explaining',
    format(sum(spectra.pca$pca.summary[selected=='Selected']$propvar),digits=2),
    'of the original sample variance')
cat('\n\tGenerating scree plot of PCA variances')
p1<-plot_elbow(spectra.pca$pca.summary)

# Save rotation matrix from PCA for selected PCs and per-gene centers
temp.dt<-dcast(final.molten[!gene_status%in%c('report') & status=='use'],sample_id~gene_id,value.var='adjlogcpkmed')
rotation.dt<-merge(data.table(mean=colMeans(temp.dt[,-1]),gene_id=names(colMeans(temp.dt[,-1]))),
      data.table(spectra.pca$pca$rotation[,1:selected.pc.count],keep.rownames=TRUE),by.x='gene_id',by.y='rn')

# Calculate scores in all samples
# Centering on original genes
temp2.dt<-dcast(final.molten[!gene_status%in%c('report')],sample_id~gene_id,value.var='adjlogcpkmed')
pc.scores<-cbind(temp2.dt[,1],scale(data.matrix(temp2.dt[,-1]),center=rotation.dt$mean,scale=FALSE)%*%data.matrix(rotation.dt[,-c(1:2)]))
pca.dt<-merge(batch.dt,pc.scores,by='sample_id')

pca_filename<-paste0(OUTPUT_DIR,'/spectra_scores-batch_variables-',Sys.Date(),'.csv')
fwrite(pca.dt,file=pca_filename)
cat('\n\tSpectra scores and batch variables written to:',pca_filename)

rotation_filename<-paste0(OUTPUT_DIR,'/spectra_rotation_matrix-gene_centers-',Sys.Date(),'.csv')
fwrite(rotation.dt,file=rotation_filename)
cat('\n\tSpectra rotation matrix and gene centers written to:',rotation_filename)

pcvar_filename<-paste0(OUTPUT_DIR,'/pca_details-',Sys.Date(),'.csv')
fwrite(spectra.pca$pca.summary,file=pcvar_filename)
cat('\n\tPCA variances and details written to:',pcvar_filename)

expr_filename<-paste0(OUTPUT_DIR,'/normalized_expression-',Sys.Date(),'.csv.gz')
fwrite(final.molten,file=expr_filename)
cat('\n\tNormalized expression data written to:',expr_filename)

scree_plot_filename=paste0(OUTPUT_DIR,'/pca_scree_plot-',Sys.Date(),'.pdf')
ggsave(scree_plot_filename,p1,height=8,width=8)
cat('\n\tSpectra scree plot written to:',scree_plot_filename)

cat('\nFinished.\n')