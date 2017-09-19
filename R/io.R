options(stringsAsFactors=F)

#' @export
read.table.smart <- function(path, ...) {
    fields <- list(...)    
    conn <- file(path)
    header <- readLines(conn, n=1)
    close(conn)
    
    sep <- "\t"
    for (s in c("\t", " ", ",")) {
        if (grepl(s, header)) {
            sep <- s
            break
        } 
    }
    
    res <- as.data.table(read.table(path, sep=sep, header=T, stringsAsFactors=F))
    
    oldnames <- character(0)
    newnames <- character(0)
    
    for (field in names(fields)) {        
        if (field %in% colnames(res)) {
            next
        }
        
        z <- na.omit(
            match(
                tolower(c(field, fields[[field]])),
                tolower(colnames(res))))
        if (length(z) == 0) {
            next
        }
        
        oldnames <- c(oldnames, colnames(res)[z])
        newnames <- c(newnames, field)
    }
        
    setnames(res, oldnames, newnames)
    res
}



#' @export
read.tsv <- function(file, header=T, sep="\t", quote="", comment.char="", check.names=FALSE, ...) {
    args <- list(...)
    res <- read.table(file, header=header, sep=sep, quote=quote, 
               comment.char=comment.char, check.names=check.names,
               stringsAsFactors=FALSE,
               ...)
    if ((!"row.names" %in% names(args)) && (colnames(res)[1] == "")) {
        rownames(res) <- res[, 1]
        res[[1]] <- NULL
    }
    res
}

#' @export
write.tsv <- function(table, dir, file=NULL, gzip=FALSE, row.names=NA, col.names=NA, ...) {
    name <- deparse(substitute(table))
    table <- as.data.frame(table) 
    
    if (is.null(file)) {
        file <- file.path(dir, paste0(name, ".tsv", if (gzip) ".gz"))        
    }

    if (is.na(row.names)) {
        row.names <- is.character(attr(table, "row.names"))
    }

    if (!row.names && is.na(col.names)) {
        col.names=T
    }
    
    for (c in colnames(table)) {
        if (is.character(table[[c]])) {
            table[[c]] <- sub("#", "", table[[c]])            
        }
    }
    
    if (gzip) {
        file <- gzfile(file, "w")
    }
    write.table(table, file, quote=F,
                row.names=row.names, col.names=col.names, sep="\t")
    if (gzip) {
        close(file)
    }
}

#' Reads gtf file to data.table
#' @param file Path to file
#' @import data.table
#' @export
read.gtf <- function(file, attrs.to.extract=c("gene_id", "transcript_id", "gene_type", "gene_name"),
                     features.to.extract=NULL) {
    res <- fread(file, header=F, col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

    if (!is.null(features.to.extract)) {
        res <- res[feature %in% features.to.extract,]
    }
                 
    attrlist <- strsplit(res$attribute, "; *")
    attrlist_length <- sapply(attrlist, length)
    attrtable <- data.table(rn=rep(seq_along(attrlist), attrlist_length), raw=unlist(attrlist))
    attrtable[, name := gsub(" .*", "", raw)]
    attrtable[, value := gsub(".* ", "", raw)]
    attrtable[, value := gsub('^"(.*)"$', "\\1", value)]
    attrtable[, raw := NULL]
    all_attrs <- unique(attrtable$name)
    attrmatrix <- matrix(nrow = length(attrlist), ncol=length(all_attrs))
    colnames(attrmatrix) <- all_attrs
    attrtable[, name := match(name, all_attrs)]
    attrmatrix[cbind(attrtable$rn, attrtable$name)] <- attrtable$value
    res[, attribute := NULL]
    
    res <- cbind(res, attrmatrix)
    res
}
