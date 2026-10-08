args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=globalenv())
cat('Loading original author workspace\n');flush.console()
objects <- load(args[1], envir=e)
cat('Objects:', paste(objects,collapse=', '),'\n');flush.console()
writeLines(objects,file.path(args[2],'author_workspace_objects.txt'))
stopifnot('snRNA_all_nuclei' %in% objects)
x <- e[['snRNA_all_nuclei']]
cat('Reference class:',class(x),'\n');flush.console()
if(isS4(x) && 'meta.data' %in% slotNames(x)) {
    m <- slot(x,'meta.data')
} else if(is.list(x) && !is.null(x$meta.data)) {
    m <- x$meta.data
} else stop('Original reference metadata structure unsupported; do not infer labels')
stopifnot(is.data.frame(m),'Cluster' %in% names(m),!anyDuplicated(rownames(m)))
cat('Metadata:',nrow(m),'rows; columns:',paste(names(m),collapse=', '),'\n');flush.console()
m$author_barcode <- rownames(m)
write.csv(m,file.path(args[2],'GSE200984_author_all_nuclei_metadata.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(args[2],'author_label_export_R_environment.txt'))
cat('Original label export completed\n');flush.console()
