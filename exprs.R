library(Biobase)
library(data.table)

# exprs.file <- "./work/data/mmp_artyomov_1195/mmu.mp.full.refseq.htcnt"
# pdata.file <- "./work/data/mmp_artyomov_1195/conditions.tsv"

makeExpressionSetFromFile <- function(
    exprs.file, 
    pdata.file,
    ignore.pattern = "^__"
    ) {
    
    require(Biobase)
    exprs <- read.table(exprs.file)
    exprs <- as.matrix(exprs[grep(ignore.pattern, rownames(exprs), invert=TRUE), ])        
    pdata <- read.table(pdata.file, header=T, row.names=1)
    pdata <- pdata[colnames(exprs), , drop=F]
    meta <- data.frame(labelDescription = colnames(pdata))
    rownames(meta) <- colnames(pdata)
    
    pdata <- new("AnnotatedDataFrame", data=pdata, varMeta=meta)
    
    eSet <- ExpressionSet(exprs, phenoData=pdata)
    eSet
}

DESeqDataSetFromExpressionSet <- function(eSet, design) {
    DESeqDataSetFromMatrix(exprs(eSet), pData(eSet), design)
}

exprsFromSpotfireDir <- function(met.data.dir) {
    files <- list.files(met.data.dir, pattern="forSpotfire.txt", full.names=T, recursive=T)
    suffixes <- sapply(files, function(f) basename(dirname(f)))
    
    stopifnot(all(!duplicated(suffixes)))
    
    raw.list <- list()
    
    for (f in files) {
        message(sprintf("Processing file %s", f))
        suffix <- suffixes[[f]]
        raw.list[[suffix]] <- exprsFromSpotfire(met.data.file=f, ion.suffix=suffix)    
    }
    
    ion.exprs <- do.call(rbind.fill, raw.list)
    rownames(ion.exprs) <- unlist(lapply(raw.list, rownames))
    
    sample.metas <- lapply(raw.list, function(x) attr(x, "sample.meta"))
    
    #     cleangroupname <- sample.metas[[1]]$cleangroupname
    #     names(cleangroupname) <- rownames(sample.metas[[1]])
    #     sapply(sample.metas, function(m) all(m$cleangroupname == cleangroupname))
    
    
    attr(ion.exprs, "ion.meta") <- do.call(rbind, unname(lapply(raw.list, function(x) attr(x, "ion.meta"))))
    attr(ion.exprs, "ion2mets") <- do.call(rbind, unname(lapply(raw.list, function(x) attr(x, "ion2mets"))))
    
    ion.exprs    
}

exprsFromSpotfire <- function(met.data.file, ion.suffix) {
    met.data <- read.csv(file=met.data.file, head=FALSE, sep="\t", stringsAsFactors=F, colClasses="character")
    met.data <- unname(t(as.matrix(met.data)))
    n <- nrow(met.data)
    m <- ncol(met.data)
    ns <- which(met.data[, 1] != "")[1]
    ms <- which(met.data[1, ] != "")[1]
    met.data[(ns+1):n, ms] <- paste0(met.data[(ns+1):n, ms], ion.suffix)
    met.exprs <- met.data[(ns+1):n, (ms+1):m]
    
    ion.meta <- as.data.frame(met.data[(ns+1):n, 1:(ms-1)], stringsAsFactors=F)
    colnames(ion.meta) <- met.data[ns, 1:(ms-1)]
    rownames(ion.meta) <- met.data[(ns+1):n, ms]
    ion.meta$mz <- as.numeric(ion.meta$mz)
    ion.meta$lowSignal <- as.logical(ion.meta$lowSignal)
    
    sample.meta <- as.data.frame(t(met.data[1:ns, (ms+1):m]), stringsAsFactors=F)
    colnames(sample.meta) <- met.data[1:ns, ms]
    rownames(sample.meta) <- paste0(sample.meta$platerow, sample.meta$platecolumn)
    
    class(met.exprs) <- "numeric"
    met.exprs <- data.frame(met.exprs)
    colnames(met.exprs) <- rownames(sample.meta)
    rownames(met.exprs) <- rownames(ion.meta)
    
    ion2mets <- ion.meta$assignedMetabolite
    names(ion2mets) <- rownames(ion.meta)
    ion2mets <- ion2mets[ion2mets != "---"]
    
    ion2mets.list <- lapply(names(ion2mets), function(ion) {
        mets.string <- ion2mets[ion]
        met.strings <- unlist(strsplit(mets.string, " /// ", fixed=T), use.names=F)
        res <- data.frame(do.call(rbind, strsplit(met.strings, " // ", fixed=T)), stringsAsFactors=F)
        colnames(res) <- c("Name", "hmdbAndAdduct", "Mass", "Score")
        res$Mass <- as.numeric(res$Mass)
        res$Score <- as.numeric(res$Score)
        hmdb.and.adduct <- do.call(rbind, strsplit(res$hmdbAndAdduct, ":", fixed=T))
        res$HMDB <- hmdb.and.adduct[,1]
        res$Adduct <- hmdb.and.adduct[,2]
        res$Ion <- ion
        res$hmdbAndAdduct <- NULL
        res
    })
    
    ion2mets.table <- do.call(rbind, ion2mets.list)
    
    #met.exprs <- merge(ion2mets.table, met.exprs)
    attr(met.exprs, "ion.meta") <- ion.meta
    attr(met.exprs, "sample.meta") <- sample.meta
    attr(met.exprs, "ion2mets") <- ion2mets.table
    met.exprs
}

rename.smart <- function(de, ...) {    
    fields <- list(...) 
    oldnames <- character(0)
    newnames <- character(0)
    
    for (field in names(fields)) {        
        if (field %in% colnames(de)) {
            next
        }
        
        z <- na.omit(
            match(
                tolower(c(field, fields[[field]])),
                tolower(colnames(de))))
        if (length(z) == 0) {
            next
        }
        z <- z[1]
        
        oldnames <- c(oldnames, colnames(de)[z])
        newnames <- c(newnames, field)
    }
    
    setnames(de, oldnames, newnames)    
}

normalizeMetDE <- function(de, org=NA, annotate=TRUE) {
    de <- as.data.table(as.data.frame(de), keep.rownames=TRUE)
    rename.smart(de, 
                 ID=c("KEGG", "HMDB", "name", "rn"),
                 pval=c("p.value", "pvalue"),
                 log2FC=c("log2foldchange", "logfc")
    )    
        
    de <- de[!duplicated(de$ID), ]
    de <- de[!is.na(de$pval), ]
    de <- as.data.table(de[order(de$pval), ])
    de
}



normalizeGeneDE <- function(de, org=NA, annotate=TRUE) {
    de <- as.data.table(as.data.frame(de), keep.rownames=TRUE)
    rename.smart(de, 
                 ID=c("gene", "entrez", "symbol", "rn"),
                 pval=c("p.value", "pvalue"),
                 log2FC=c("log2foldchange", "logfc")
    )    
    # :ToDo: it's a hack
    
    if (annotate) {
        de <- merge(as.data.frame(de), unique(reflink[, list(Entrez, symbol, product)]), by.x="ID", by.y="Entrez", all.x=T)
    }
    de <- de[!duplicated(de$ID), ]
    de <- de[!is.na(de$pval), ]
    de <- as.data.table(de[order(de$pval), ])
    de
}


read.gct <- function(gct, ...) { 
    meta <- readLines(gct, n=2)
    version <- meta[1]
    size <- as.numeric(unlist(strsplit(meta[2], "\t")))
    
    stopifnot(grepl("^#1.2", version))
    
    t <- read.tsv(gct, skip=2, nrows=size[1], row.names=1, ...)
    
    exp <- as.matrix(t[, rev(ncol(t) + 1 - seq_len(size[2]))])
    rownames(exp) <- rownames(t)
    
    fdata <- t[,seq_len(ncol(t) - size[2]), drop=F]
    fmeta <- data.frame(labelDescription = colnames(fdata))
    rownames(fmeta) <- colnames(fdata)
    
    fdata <- new("AnnotatedDataFrame", data=fdata, varMeta=fmeta)
    
    ExpressionSet(exp, featureData=fdata)
}

makeAnnotated <- function(data) {
    meta <- data.frame(labelDescription = colnames(data))
    rownames(meta) <- colnames(data)
    
    new("AnnotatedDataFrame", data=data, varMeta=meta)    
}

zScore <- function(x) {
    x.means <- apply(x, 1,mean)
    x.sds<- apply(x, 1, sd)
    res <- sweep(sweep(x, 1, x.means), 1, x.sds, "/")
    return(res)
}

normalize.rows <- function(x) { 
    x <- sweep(x, 1, apply(x, 1, min)) 
    sweep(x, 1, apply(x, 1, max), "/") 
} 

