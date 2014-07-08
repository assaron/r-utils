# exprs.file <- "./work/data/mmp_artyomov_1195/mmu.mp.full.refseq.htcnt"
# pdata.file <- "./work/data/mmp_artyomov_1195/conditions.tsv"

makeExpressionSetFromFile <- function(
    exprs.file, 
    pdata.file
    ) {
    
    require(Biobase)
    exprs <- read.table(exprs.file)
    exprs <- as.matrix(exprs)
    pdata <- read.table(pdata.file, header=T, row.names=1)
    pdata <- pdata[colnames(exprs), , drop=F]
    meta <- data.frame(labelDescription = colnames(pdata))
    rownames(meta) <- colnames(pdata)
    
    pdata <- new("AnnotatedDataFrame", data=pdata, varMeta=meta)
    
    eSet <- ExpressionSet(exprs, phenoData=pdata)
    eSet
}

