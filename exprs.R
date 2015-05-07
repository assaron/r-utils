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
    exprs <- read.tsv(exprs.file, header=T, row.names=1)
    exprs <- as.matrix(exprs[grep(ignore.pattern, rownames(exprs), invert=TRUE), ])        
    pdata <- read.tsv(pdata.file, header=T, row.names=1)
    missing <- setdiff(rownames(pdata), colnames(exprs))
    if (length(missing) > 0) {
        warning(paste0("There are missing samples: ", paste(missing, collapse=" ")))        
        
    }
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

makeAnnotated <- function(data) {
    meta <- data.frame(labelDescription = colnames(data))
    rownames(meta) <- colnames(data)
    
    new("AnnotatedDataFrame", data=data, varMeta=meta)    
}

read.gct <- function(gct, ...) { 
    meta <- readLines(gct, n=3)
    version <- meta[1]
    size <- as.numeric(unlist(strsplit(meta[2], "\t")))
    
    if (grepl("^#1.3", version)) {
        ann.col <- size[4]
        ann.row <- size[3]
    } else if (grepl("^#1.2", version)) {
        ann.col <- 1
        ann.row <- 0
    } else {
        stop("Unsupported version of gct: use 1.2 or 1.3")
    }
    
    
    
    t <- read.tsv(gct, skip=2 + 1 + ann.col, nrows=size[1], col.names=unlist(strsplit(meta[3], "\t")), row.names=1, header=F,  ...)    
    
        
    exp <- as.matrix(t[, (ann.row+1):ncol(t)])
        
    fdata <- makeAnnotated(t[,seq_len(ann.row), drop=F])

    
    if (ann.row > 0) {
        pdata.raw <- t(read.tsv(gct, skip=2+1, nrows=ann.row, header=F))
        pdata <- data.frame(pdata.raw[seq_len(ncol(exp))+1+ann.row, , drop=F])
        colnames(pdata) <- pdata.raw[1,]
        rownames(pdata) <- colnames(exp)
        pdata <- makeAnnotated(pdata)       
        
        res <- ExpressionSet(exp, featureData=fdata, phenoData=pdata)
    } else {        
        res <- ExpressionSet(exp, featureData=fdata)
    }
    
    res    
}

write.gct <- function(es, file, gzip=FALSE) {
    if (gzip) {
        con <- gzfile(file)
    } else {
        con <- file(file)
    }
    open(con, open="w")
    writeLines("#1.3", con)
    ann.col <- ncol(pData(es))
    ann.row <- ncol(fData(es))
    writeLines(sprintf("%s\t%s\t%s\t%s", nrow(es), ncol(es), ann.row, ann.col), con)
    writeLines(paste0(c("ID", colnames(fData(es)), colnames(es)), collapse="\t"), con)
    
    ann.col.table <- t(as.matrix(pData(es)))    
    ann.col.table <- cbind(matrix(rep(NA, ann.row*ann.col), nrow=ann.col), ann.col.table)
    write.table(ann.col.table, file=con, quote=F, sep="\t", row.names=T, col.names=F)                          
    write.table(cbind(fData(es), exprs(es)), file=con, quote=F, sep="\t", row.names=T, col.names=F)                          
    close(con)    
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

