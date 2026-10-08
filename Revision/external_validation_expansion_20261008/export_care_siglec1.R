suppressPackageStartupMessages(library(Matrix))
out <- file.path(Sys.getenv('JOB34_ROOT',unset='.'),'Revision/external_validation_expansion_20261008')
src <- Sys.getenv('JOB34_CARE_METADATA',unset='external_inputs/CARE')
dir.create(file.path(out,'care_metrics'),showWarnings=FALSE)
meta <- read.csv(file.path(src,'celltype_meta_data_2025_01_08.csv'))
clinical <- read.csv(file.path(src,'CARE_original_clinical.csv'),check.names=FALSE)
names(clinical) <- trimws(names(clinical))
for(k in names(clinical)) if(is.character(clinical[[k]])) clinical[[k]] <- trimws(clinical[[k]])
ages <- tapply(clinical[['Age at initial diagnosis']],clinical[['Patient ID']],function(z){z<-z[!is.na(z)];if(length(z))z[1] else NA_real_})
clinical$age_initial <- as.numeric(ages[clinical[['Patient ID']]])
clinical$eligible <- !is.na(clinical$age_initial) & clinical$age_initial>=18 & clinical[['IDH mutation status']]=='IDHwt'
clinical$condition <- ifelse(clinical[['Primary or recurrent']]=='Primary','primary',ifelse(grepl('Recurrent',clinical[['Primary or recurrent']]),'recurrent',NA))
write.csv(clinical,file.path(out,'CARE_clinical_eligibility.csv'),row.names=FALSE)
ids <- sort(intersect(unique(meta$ID),clinical[['Sample ID']][clinical$eligible & !is.na(clinical$condition)]))
genes <- c('SIGLEC1','CD68','AIF1','LST1','TYROBP','CSF1R')
programs <- list(
 Interferon_response=c('ISG15','IFIT1','IFIT2','IFIT3','IFI6','IFI27','IFI44L','MX1','OAS1','STAT1','IRF7'),
 Phagolysosome=c('LAMP1','CTSD','CTSB','FCGR1A','TYROBP','AIF1','LILRB1','CTSS'),
 Lipid_and_complement=c('APOE','LPL','TREM2','C1QA','C1QB','C1QC','GPNMB','FABP5'),
 Antigen_presentation=c('HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','B2M'),
 Immunoregulatory=c('IL10','TGFB1','CD274','ARG1','MRC1','MSR1','LGALS3'),
 Immediate_inflammatory=c('CCL3','CCL4','IL1B','TNF','NFKBIA','CXCL8','PTGS2','FOS','JUN'))
stopifnot(!anyDuplicated(meta$CellID),!anyDuplicated(clinical[['Sample ID']]))
for(id in ids) {
 dest<-file.path(out,'care_metrics',paste0(id,'.csv'))
 if(file.exists(dest)){cat(id,'cached\n');next}
 paths<-list.files(Sys.getenv('JOB34_CARE_RAW',unset='external_inputs/CARE/raw'),pattern=paste0('_',id,'_umi_counts.RDS.gz$'),full.names=TRUE)
 if(length(paths)!=1)stop(paste('missing/ambiguous',id))
 x<-readRDS(gzcon(gzfile(paths[1],open='rb')))
 if(!inherits(x,'sparseMatrix'))x<-Matrix(x,sparse=TRUE)
 stopifnot(inherits(x,'sparseMatrix'),!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(x@x>=0),all(genes%in%rownames(x)))
 m<-meta[meta$ID==id,];j<-match(m$CellID,colnames(x));stopifnot(!anyNA(j))
 x<-x[,j,drop=FALSE];lib<-colSums(x);ng<-colSums(x>0);mt<-colSums(x[grepl('^MT-',rownames(x)),,drop=FALSE])/pmax(lib,1)
 keep<-lib>=500 & ng>=300 & ng<=8000 & mt<=.2
 rows<-list()
 for(group in c('TAM','All_nuclei')) {
  ix<-which(keep & (group=='All_nuclei' | m$CellType=='TAM'))
  if(!length(ix))next
  y<-log1p(1e4*as.matrix(x[genes,ix,drop=FALSE])/rep(lib[ix],each=length(genes)))
  d<-data.frame(sample=id,group=group,n_cells=length(ix),n_author=nrow(m),total_QC=sum(keep),myeloid_fraction=sum(keep&m$CellType=='TAM')/sum(keep),library_counts=sum(lib[ix]),SIGLEC1_detection=mean(y['SIGLEC1',]>0),mean_ln_UMI=mean(log1p(lib[ix])),mean_ln_genes=mean(log1p(ng[ix])))
  for(g in genes){d[[g]]<-mean(y[g,]);d[[paste0('counts_',g)]]<-sum(x[g,ix])}
  for(prog in names(programs)) {
   gs<-intersect(programs[[prog]],rownames(x));d[[paste0('coverage_',prog)]]<-length(gs)/length(programs[[prog]])
   if(length(gs)<ceiling(.8*length(programs[[prog]])))stop(paste('low gene coverage',id,prog))
   score<-colMeans(log1p(1e4*as.matrix(x[gs,ix,drop=FALSE])/rep(lib[ix],each=length(gs))))
   d[[prog]]<-mean(score)
   r<-if(sd(y['SIGLEC1',])>0 && sd(score)>0)cor(y['SIGLEC1',],score,method='spearman') else NA_real_
   d[[paste0('rho_',prog)]]<-r
   if(length(ix)>5 && sd(y['SIGLEC1',])>0 && sd(score)>0) {
    design<-cbind(1,rank(log1p(lib[ix])))
    a<-lm.fit(design,rank(y['SIGLEC1',]))$residuals;b<-lm.fit(design,rank(score))$residuals
    d[[paste0('depth_rho_',prog)]]<-if(sd(a)>0&&sd(b)>0)cor(a,b) else NA_real_
   } else d[[paste0('depth_rho_',prog)]]<-NA_real_
  }
  rows[[length(rows)+1]]<-d
 }
 write.csv(do.call(rbind,rows),dest,row.names=FALSE)
 cat(id,'author',nrow(m),'QC',sum(keep),'TAM',sum(keep&m$CellType=='TAM'),'\n');flush.console()
 rm(x,y);gc(verbose=FALSE)
}
files<-list.files(file.path(out,'care_metrics'),pattern='csv$',full.names=TRUE)
write.csv(do.call(rbind,lapply(files,read.csv)),file.path(out,'CARE_specimen_metrics.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'CARE_R_environment.txt'))
cat('CARE_SIGLEC1_EXPORT_COMPLETE\n')
